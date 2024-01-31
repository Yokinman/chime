use proc_macro::*;
use quote::*;
use syn::parse::ParseStream;
use syn::Token;

struct FluxParse {
	kind_type: syn::Type,
	change_expr: syn::Expr,
	crate_path: syn::Path,
}

impl syn::parse::Parse for FluxParse {
	fn parse(input: ParseStream) -> syn::Result<Self> {
		 // Parse Main Change Expression:
		let kind_type = syn::Type::parse(input)?;
		input.parse::<Token![=]>()?;
		let change_expr = syn::Expr::parse(input)?;
		
		 // Parse Crate Name Override:
		let mut crate_path: syn::Path = syn::parse_quote!{::chime};
		while !input.is_empty() {
			input.parse::<Token![,]>()?;
			let meta = syn::MetaNameValue::parse(input)?;
			if meta.path.is_ident("crate") {
				if let syn::Expr::Lit(syn::ExprLit { attrs, lit: syn::Lit::Str(s) }) = meta.value {
					if !attrs.is_empty() {
						panic!("no attributes allowed for crate path");
					}
					crate_path = s.parse()?;
				} else {
					panic!(
						"unexpected meta expression -> {}",
						meta.value.to_token_stream()
					)
				}
			} else {
				panic!(
					"invalid meta name -> {}",
					meta.path.to_token_stream()
				);
			}
		}
		
		Ok(FluxParse {
			kind_type,
			change_expr,
			crate_path,
		})
	}
}

fn contextualize(
	expr: &mut syn::Expr,
	value_ident: &mut Option<syn::Ident>,
	fields: &syn::Fields,
	flux: &proc_macro2::TokenStream,
) {
	match expr {
		syn::Expr::Paren(syn::ExprParen { expr, .. }) => {
			contextualize(expr, value_ident, fields, flux);
		},
		// syn::Expr::Match(syn::ExprMatch { expr, arms, .. }) => {
		// 	contextualize(expr, value_ident, fields, flux);
		// 	for arm in arms {
		// 		contextualize(arm.body.as_mut(), value_ident, fields, flux);
		// 	}
		// },
		syn::Expr::MethodCall(syn::ExprMethodCall { receiver, /*args,*/ .. }) => {
			contextualize(receiver, value_ident, fields, flux);
			// for arg in args {
			// 	contextualize(arg, value_ident, fields, flux);
			// }
		},
		syn::Expr::Binary(syn::ExprBinary { left, right, .. }) => {
			contextualize(left,  value_ident, fields, flux);
			contextualize(right, value_ident, fields, flux);
		},
		syn::Expr::Unary(syn::ExprUnary { expr, .. }) => {
			contextualize(expr, value_ident, fields, flux);
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
				contextualize(&mut path, &mut None, /*&mut Default::default(),*/ fields, flux);
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
			if let Some(_) = fields.iter().find(|field| path.is_ident(
				&field.ident.as_ref().unwrap().to_string()
			)) {
				let ident = path.get_ident();
				let mut expr_field: syn::ExprField = syn::parse_quote!{
					self.#ident
				};
				expr_field.attrs = std::mem::take(attrs);
				*expr = expr_field.into();
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
	let (kind_type, value_ident, change_expr, flux) = match syn::parse(arg_stream) {
		Ok(FluxParse { kind_type, mut change_expr, crate_path: flux }) => {
			let flux = flux.to_token_stream();
			let mut value_ident = None;
			
			contextualize(&mut change_expr, &mut value_ident, &item.fields, &flux);
			
			let value_ident = if let Some(value_expr) = value_ident {
				value_expr
			} else {
				panic!("no value identifier found (wrap one in curly brackets)");
			};
			
			(kind_type, value_ident, change_expr, flux)
		}
		Err(err) => panic!("{}", err),
	};
	
	// println!("item        : {}", item.to_token_stream().to_string());
	// println!("kind_type   : {}", kind_type.to_token_stream().to_string());
	// println!("value_ident : {}", value_ident.to_token_stream().to_string());
	// println!("change_expr : {}", change_expr.to_token_stream().to_string());
	
	 // Generate `Flux` Variant of Struct:
	let mut flux_item = item.clone();
	let flux_ident = syn::Ident::new(
		format!("{}Flux", item.ident).as_str(),
		proc_macro2::Span::mixed_site()
	);
	flux_item.ident = flux_ident.clone();
	let mut moment_fields = Default::default();
	let mut flux_fields = Default::default();
	for field in &mut flux_item.fields {
		let ident = field.ident.as_ref();
		if ident == Some(&value_ident) {
			let field_ty = std::mem::replace(&mut field.ty, syn::parse_quote!{
				#flux::Constant<<<Self as #flux::Flux>::Kind
					as #flux::kind::FluxKind>::Value>
			});
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::linear::LinearIso::<#field_ty>
					::map(self.value(time)),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: #flux::Moment::to_flux(
					&#flux::linear::LinearIsoInv
						::<<<Self::Flux as #flux::Flux>::Kind
							as #flux::kind::FluxKind>::Value>
						::inv_map(self.#ident),
					time
				),
			};
		} else {
			let field_ty = &field.ty;
			field.ty = syn::parse_quote!{
				<#field_ty as #flux::Moment>::Flux
			};
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::Flux::to_moment(&self.#ident, time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: #flux::Moment::to_flux(&self.#ident, time),
			};
		}
	}
	
	 // Generate Implementation:
	let ident = item.ident.clone();
	let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();
	let flux_value_impl = quote::quote!{
		#item
		
		impl #impl_generics #flux::Moment for #ident #ty_generics #where_clause {
			type Flux = #flux_ident #ty_generics;
			fn to_flux(&self, time: #flux::time::Time) -> Self::Flux {
				#flux_ident { #flux_fields }
			}
		}
		
		// ??? Unsure if `item.attrs` and `item.vis` should be inherited by
		// `flux_item` as they currently are.
		// ??? Unsure if `flux_item` should be placed in a private module, so
		// it's only accessible through `flux::Moment::Flux`.
		#flux_item
		
		impl #impl_generics #flux::Flux for #flux_ident #ty_generics #where_clause {
			type Moment = #ident #ty_generics;
			type Kind = #kind_type;
			fn base_value(&self) -> <Self::Kind as #flux::kind::FluxKind>::Value {
				self.#value_ident.base_value()
			}
			fn base_time(&self) -> #flux::time::Time {
				self.#value_ident.base_time()
			}
			fn change<'a>(&self, accum: <Self::Kind as #flux::kind::FluxKind>::Accum<'a>)
				-> <Self::Kind as #flux::kind::FluxKind>::OutAccum<'a>
			{
				#change_expr
			}
			fn to_moment(&self, time: #flux::time::Time) -> Self::Moment {
				#ident { #moment_fields }
			}
		}
	};
	
	// println!("flux_value_impl: {}", flux_value_impl.to_string());
	
	flux_value_impl.into()
}