use proc_macro::*;
use quote::*;
use syn::parse::ParseStream;
use syn::Token;

struct FluxParse {
	kind_type: syn::Type,
	change_expr: syn::Expr,
}

impl syn::parse::Parse for FluxParse {
	fn parse(input: ParseStream) -> syn::Result<Self> {
		let kind_type = syn::Type::parse(input)?;
		input.parse::<Token![=]>()?;
		let change_expr = syn::Expr::parse(input)?;
		Ok(FluxParse {
			kind_type,
			change_expr,
		})
	}
}

fn contextualize(
	expr: &mut syn::Expr,
	value_ident: &mut Option<syn::Ident>,
	out_accum: &mut proc_macro2::TokenStream,
	fields: &syn::Fields,
) {
	match expr {
		syn::Expr::Paren(syn::ExprParen { expr, .. }) => {
			contextualize(expr, value_ident, out_accum, fields);
		},
		
		 // Log Operation Outputs:
		syn::Expr::Binary(syn::ExprBinary { left, op, right, .. }) => {
			let mut lhs = out_accum.clone();
			let mut rhs = std::mem::take(out_accum);
			contextualize(left,  value_ident, &mut lhs, fields);
			contextualize(right, value_ident, &mut rhs, fields);
			let op_trait = match op {
				syn::BinOp::Add(_)    => quote::quote!{Add},
				syn::BinOp::Sub(_)    => quote::quote!{Sub},
				syn::BinOp::Mul(_)    => quote::quote!{Mul},
				syn::BinOp::Div(_)    => quote::quote!{Div},
				syn::BinOp::Rem(_)    => quote::quote!{Rem},
				syn::BinOp::BitAnd(_) => quote::quote!{BitAnd},
				syn::BinOp::BitOr(_)  => quote::quote!{BitOr},
				syn::BinOp::BitXor(_) => quote::quote!{BitXor},
				syn::BinOp::Shl(_)    => quote::quote!{Shl},
				syn::BinOp::Shr(_)    => quote::quote!{Shr},
				e => panic!(
					"unexpected binary expression -> {}",
					e.to_token_stream()
				)
			};
			*out_accum = quote::quote!{
				<#lhs as std::ops::#op_trait<#rhs>>::Output
			};
		},
		syn::Expr::Unary(syn::ExprUnary { op, expr, .. }) => {
			let mut value = out_accum.clone();
			contextualize(expr, value_ident, &mut value, fields);
			let op_trait = match op {
				syn::UnOp::Neg(_) => quote::quote!{Neg},
				syn::UnOp::Not(_) => quote::quote!{Not},
				e => panic!(
					"unexpected unary expression -> {}",
					e.to_token_stream()
				)
			};
			*out_accum = quote::quote!{
				<#value as std::ops::#op_trait>::Output
			};
		},
		
		 // Disambiguate `Per::per` Method Call:
		syn::Expr::MethodCall(syn::ExprMethodCall {
			attrs, receiver, method, turbofish: None, args, ..
		}) => {
			if *method != "per" || args.len() != 1 {
				panic!("unexpected method call (can only call `Per::per`)");
			}
			let mut ty = std::mem::take(out_accum);
			contextualize(receiver, value_ident, &mut ty, fields);
			let mut call: syn::ExprCall = syn::parse_quote!{
				flux::Per::#method(&#receiver, #args)
			};
			call.attrs = std::mem::take(attrs);
			*expr = syn::Expr::Call(call);
			*out_accum = quote::quote!(flux::Change<'a, <#ty as flux::Moment>::Flux>);
		},
		
		 // Retrieve Value Expression (`{value}`):
		syn::Expr::Block(syn::ExprBlock { block, attrs, label: None }) => {
			if block.stmts.len() != 1 {
				panic!("expected a value identifier");
			}
			if let Some(syn::Stmt::Expr(syn::Expr::Path(e), None)) = block.stmts.last_mut() {
				 // Store Value Expression:
				if let Some(value_expr) = value_ident {
					panic!(
						"found multiple value expressions\n1) {}\n2) {}",
						value_expr.to_token_stream(),
						e.to_token_stream(),
					);
				}
				let mut path = syn::Expr::Path(e.clone());
				contextualize(&mut path, &mut None, &mut Default::default(), fields);
				*value_ident = match path {
					syn::Expr::Field(syn::ExprField {
						member: syn::Member::Named(ident),
						..
					}) => Some(ident),
					_ => unreachable!()
				};
				
				 // Replace w/ Identifier:
				let mut accum: syn::ExprPath = syn::parse_str("accum").unwrap();
				accum.attrs = std::mem::take(attrs);
				*expr = accum.into();
			} else {
				panic!(
					"unexpected block expression -> {}",
					block.to_token_stream()
				)
			}
		},
		
		 // Convert Identifiers to `self` Fields:
		syn::Expr::Path(syn::ExprPath { attrs, path, .. }) => {
			if path.is_ident("self") {
				panic!("can't use self here");
			}
			if let Some(field) = fields.iter().find(|field| path.is_ident(
				&field.ident.as_ref().unwrap().to_string()
			)) {
				let ident = path.get_ident();
				let mut expr_field: syn::ExprField = syn::parse_quote!{
					self.#ident
				};
				expr_field.attrs = std::mem::take(attrs);
				*expr = expr_field.into();
				*out_accum = field.ty.to_token_stream();
			} else {
				panic!("can only refer to the type's fields");
			}
		},
		
		e => panic!(
			"unexpected expression -> {}",
			e.to_token_stream()
		)
	}
}

/// - Args: `KindType = ChangeExpr`.
/// - KindType: Any type that implements `FluxKind`.
/// - ChangeExpr: Any expression made up of only binary/unary operations,
///   `per` method calls, parenthesized expressions, and identifiers.
///   - The identifiers can only be names of the struct's fields.
///   - Can also have block expressions, but only ones that contain a single
///     identifier. This signifies the core "value" field of the struct.
#[proc_macro_attribute]
pub fn flux(arg_stream: TokenStream, item_stream: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_stream) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	 // Parse Attribute Arguments:
	let (kind_type, value_ident, change_expr, out_accum) = match syn::parse(arg_stream) {
		Ok(FluxParse { kind_type, mut change_expr }) => {
			let mut value_ident = None;
			let mut out_accum = quote::quote!{
				<Self::Kind as flux::FluxKind>::Accum<'a>
			};
			
			contextualize(&mut change_expr, &mut value_ident, &mut out_accum, &item.fields);
			
			let value_ident = if let Some(value_expr) = value_ident {
				value_expr
			} else {
				panic!("no value identifier found (wrap one in curly brackets)");
			};
			
			(kind_type, value_ident, change_expr, out_accum)
		},
		Err(err) => panic!("{}", err),
	};
	
	// println!("item        : {}", item.to_token_stream().to_string());
	// println!("kind_type   : {}", kind_type.to_token_stream().to_string());
	// println!("value_ident : {}", value_ident.to_token_stream().to_string());
	// println!("change_expr : {}", change_expr.to_token_stream().to_string());
	// println!("out_accum   : {}", out_accum);
	
	 // Generate `Flux` Variant of Struct:
	let mut flux_item = item.clone();
	let flux_ident = syn::Ident::new(
		format!("{}Flux", item.ident).as_str(),
		item.ident.span() // I don't understand span enough to know if this makes sense
	);
	flux_item.ident = flux_ident.clone();
	let mut moment_fields = Default::default();
	let mut flux_fields = Default::default();
	for field in &mut flux_item.fields {
		let field_ty = &field.ty;
		let ident = field.ident.as_ref();
		if ident == Some(&value_ident) {
			field.ty = syn::parse_quote!{
				flux::FluxValue<<<Self as flux::Flux>::Kind as flux::FluxKind>::Value>
			};
			moment_fields = quote::quote!{
				#moment_fields
				#ident: flux::linear::LinearIso::map(self.value_at(time)),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: flux::FluxValue::new(
					time,
					flux::linear::InvLinearIso::inv_map(self.#ident)
				),
			};
		} else {
			field.ty = syn::parse_quote!{
				<#field_ty as flux::Moment>::Flux
			};
			moment_fields = quote::quote!{
				#moment_fields
				#ident: flux::Flux::at(&self.#ident, time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: flux::Moment::to_flux(self.#ident, time),
			};
		}
	}
	
	 // Generate Implementation:
	let ident = item.ident.clone();
	let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();
	let flux_value_impl = quote::quote!{
		#item
		
		impl #impl_generics flux::Moment for #ident #ty_generics #where_clause {
			type Flux = #flux_ident #ty_generics;
			fn to_flux(self, time: flux::Time) -> Self::Flux {
				#flux_ident { #flux_fields }
			}
		}
		
		// ??? Unsure if `item.attrs` should be inherited by `flux_item`, as it
		// currently does.
		// ??? Unsure if `flux_item` should be placed in a private module, so
		// it's only accessible through `flux::Moment::Flux`.
		#flux_item
		
		impl #impl_generics flux::Flux for #flux_ident #ty_generics #where_clause {
			type Moment = #ident #ty_generics;
			type Kind = #kind_type;
			type OutAccum<'a> = #out_accum where <Self::Kind as flux::FluxKind>::Value: 'a;
			fn value(&self) -> <Self::Kind as flux::FluxKind>::Value {
				*self.#value_ident
			}
			fn time(&self) -> Time {
				self.#value_ident.time()
			}
			fn change<'a>(&self, accum: <Self::Kind as flux::FluxKind>::Accum<'a>) -> Self::OutAccum<'a> {
				#change_expr
			}
			fn at(&self, time: flux::Time) -> Self::Moment {
				#ident { #moment_fields }
			}
		}
	};
	
	// println!("flux_value_impl: {}", flux_value_impl.to_string());
	
	flux_value_impl.into()
}