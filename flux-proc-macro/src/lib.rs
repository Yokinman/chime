use proc_macro::*;
use std::collections::HashSet;
use quote::*;
use syn::ext::IdentExt;
use syn::parse::ParseStream;
use syn::Token;

struct FluxParse {
	kind_type: Option<syn::Type>,
	value_expr: Option<syn::Expr>,
	change_expr: Option<syn::Expr>,
	flux_type: Option<syn::Type>,
	crate_path: Option<syn::Path>,
}

impl syn::parse::Parse for FluxParse {
	fn parse(input: ParseStream) -> syn::Result<Self> {
		let mut kind_type = None;
		let mut value_expr = None;
		let mut change_expr = None;
		let mut flux_type = None;
		let mut crate_path = None;
		
		let mut prev_paths = HashSet::new();
		while !input.is_empty() {
			let field = syn::Path::parse(input)?;
			let ident = field.require_ident()?;
			input.parse::<Token![=]>()?;
			
			 // No Repeats:
			if prev_paths.contains(ident) {
				panic!("must only define each field once");
			}
			prev_paths.insert(ident.clone());
			
			match ident.to_string().as_str() {
				"kind" => {
					kind_type = Some(syn::Type::parse(input)?);
				},
				"value" => {
					let ident = syn::Ident::parse(input)?;
					if ident == "self" {
						panic!("expected value to be an identifier for a field");
					}
					value_expr = Some(syn::parse_quote!(#ident));
				},
				"change" => {
					let value = syn::Expr::parse(input)?;
					if let syn::Expr::Closure(_) = value {
						change_expr = Some(value);
					} else {
						panic!("change expression must be a closure");
					}
				},
				"flux" => {
					flux_type = Some(syn::Type::parse(input)?);
				},
				"crate" => {
					crate_path = Some(syn::Path::parse(input)?);
				},
				ident => panic!("invalid identifier -> {}", ident)
			}
			
			if input.is_empty() {
				break
			}
			
			input.parse::<Token![,]>()?;
		}
		
		Ok(FluxParse {
			kind_type,
			value_expr,
			change_expr,
			flux_type,
			crate_path,
		})
	}
}

fn fetch_idents(
	expr: &syn::Expr,
	idents: &mut HashSet<proc_macro2::Ident>,
	cond_idents: &mut HashSet<proc_macro2::Ident>,
) {
	match expr {
		syn::Expr::Paren(syn::ExprParen { expr, .. }) => {
			fetch_idents(expr, idents, cond_idents);
		},
		syn::Expr::Call(syn::ExprCall { func, args, .. }) => {
			let mut sub_cond_idents = HashSet::new();
			fetch_idents(func, cond_idents, &mut sub_cond_idents);
			*cond_idents = cond_idents.union(&sub_cond_idents).map(|x| x.clone()).collect();
			for arg in args {
				fetch_idents(arg, idents, cond_idents);
			}
		},
		syn::Expr::MethodCall(syn::ExprMethodCall { receiver, args, .. }) => {
			fetch_idents(receiver, idents, cond_idents);
			for arg in args {
				fetch_idents(arg, idents, cond_idents);
			}
		},
		syn::Expr::Field(syn::ExprField { base, .. }) => {
			fetch_idents(base, idents, cond_idents);
		},
		syn::Expr::Lit(syn::ExprLit { .. }) => {},
		syn::Expr::Array(syn::ExprArray { elems, .. }) => {
			for elem in elems {
				fetch_idents(elem, idents, cond_idents);
			}
		},
		syn::Expr::Closure(syn::ExprClosure { body, .. }) => {
			fetch_idents(body, idents, cond_idents);
		},
		syn::Expr::Repeat(syn::ExprRepeat { expr, len, .. }) => {
			let mut sub_cond_idents = HashSet::new();
			fetch_idents(len, cond_idents, &mut sub_cond_idents);
			*cond_idents = cond_idents.union(&sub_cond_idents).map(|x| x.clone()).collect();
			fetch_idents(expr, idents, cond_idents);
		},
		syn::Expr::Struct(syn::ExprStruct { fields, rest, .. }) => {
			for field in fields {
				fetch_idents(&field.expr, idents, cond_idents);
			}
			if let Some(rest) = rest {
				fetch_idents(rest, idents, cond_idents);
			}
		},
		syn::Expr::Tuple(syn::ExprTuple { elems, .. }) => {
			for elem in elems {
				fetch_idents(elem, idents, cond_idents);
			}
		},
		syn::Expr::Binary(syn::ExprBinary { left, right, .. }) => {
			fetch_idents(left,  idents, cond_idents);
			fetch_idents(right, idents, cond_idents);
		},
		syn::Expr::Index(syn::ExprIndex { expr, index, .. }) => {
			let mut sub_cond_idents = HashSet::new();
			fetch_idents(index, cond_idents, &mut sub_cond_idents);
			*cond_idents = cond_idents.union(&sub_cond_idents).map(|x| x.clone()).collect();
			fetch_idents(expr, idents, cond_idents);
		},
		syn::Expr::Range(syn::ExprRange { .. }) => {},
		syn::Expr::Reference(syn::ExprReference { expr, .. }) |
		syn::Expr::Unary(syn::ExprUnary { expr, .. }) |
		syn::Expr::Cast(syn::ExprCast { expr, .. }) |
		syn::Expr::Try(syn::ExprTry { expr, .. }) => {
			fetch_idents(expr, idents, cond_idents);
		},
		syn::Expr::If(syn::ExprIf { cond, then_branch, else_branch, .. }) => {
			let mut sub_cond_idents = HashSet::new();
			fetch_idents(cond, cond_idents, &mut sub_cond_idents);
			*cond_idents = cond_idents.union(&sub_cond_idents).map(|x| x.clone()).collect();
			for stmt in &then_branch.stmts {
				match stmt {
					syn::Stmt::Expr(ref expr, None) => {
						fetch_idents(expr, idents, cond_idents);
					},
					s => panic!(
						"must only use expressions -> {}",
						s.to_token_stream()
					),
				}
			}
			if let Some((_, expr)) = else_branch {
				fetch_idents(expr, idents, cond_idents);
			}
		},
		syn::Expr::Match(syn::ExprMatch { expr, arms, .. }) => {
			let mut sub_cond_idents = HashSet::new();
			fetch_idents(expr, cond_idents, &mut sub_cond_idents);
			*cond_idents = cond_idents.union(&sub_cond_idents).map(|x| x.clone()).collect();
			for arm in arms {
				fetch_idents(&arm.body, idents, cond_idents);
			}
		},
		syn::Expr::Block(syn::ExprBlock { block, .. }) => {
			for stmt in &block.stmts {
				match stmt {
					syn::Stmt::Expr(ref expr, None) => {
						fetch_idents(expr, idents, cond_idents);
					},
					s => panic!(
						"must only use expressions -> {}",
						s.to_token_stream()
					),
				}
			}
		},
		
		 // Fetch Identifier:
		syn::Expr::Path(syn::ExprPath { path, .. }) => {
			if path.segments.len() == 1 && !path.is_ident("self") {
				idents.insert(path.get_ident()
					.expect("identifier must exist")
					.unraw());
			}
		},
		
		e => panic!(
			"unexpected expression -> {}",
			e.to_token_stream()
		)
	}
}

/// Args: `kind_type="Type"`, `value=expr`, `change=closure`, `crate="path"`.
/// - `kind_type`: Any type that implements `FluxKind`.
/// - `value_expr`: An expression that describes the initial value.
/// - `change_expr`: A closure that applies change over time to an accumulator
///   and returns the result. Defaults to an identity function.
/// - `crate`: The path to the `chime` crate. Defaults to `"::chime"`.
#[proc_macro_attribute]
pub fn flux(arg_stream: TokenStream, item_stream: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_stream) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	 // Parse Attribute Arguments:
	let (
		kind_type,
		value_block,
		value_idents,
		change_block,
		change_idents,
		flux_type,
		flux,
	) = match syn::parse(arg_stream) {
		Ok(FluxParse { kind_type, value_expr, change_expr, flux_type, crate_path }) => {
			let kind_type = kind_type.expect("must specify kind type (kind_type = Type)");
			let value_expr = value_expr.expect("must specify value expression (value = expr)");
			let change_expr = change_expr.unwrap_or_else(|| syn::parse_quote!{|x| x});
			let flux_type = flux_type.unwrap_or_else(|| {
				let ident = item.ident.clone();
				let generics = item.generics.clone();
				syn::parse_quote!{#ident #generics}
			});
			let crate_path = crate_path.unwrap_or_else(|| syn::parse_quote!{::chime});
			
			let mut value_idents = HashSet::new();
			let mut value_cond_idents = HashSet::new();
			let mut change_idents = HashSet::new();
			let mut change_cond_idents = HashSet::new();
			
			// !!! Assert no overlap between `value_idents` and `change_idents`.
			
			fetch_idents(&value_expr, &mut value_idents, &mut value_cond_idents);
			assert_eq!(value_idents.len(), 1, "must specify a single value identifier without `self.`");
			fetch_idents(&change_expr, &mut change_idents, &mut change_cond_idents);
			
			let mut value_block: syn::Block = syn::parse_quote!{{#value_expr}};
			let mut change_block: syn::Block = syn::parse_quote!{{
				(#change_expr)(kind)
			}};
			
			 // Convenient Identifiers:
			let has_ident = |ident: &&syn::Ident| {
				item.fields.iter().any(|field|
					**ident == field.ident.as_ref()
						.expect("identifier must exist")
						.unraw()
				)
			};
			for ident in change_idents.iter().filter(has_ident) {
				change_block.stmts.insert(0, syn::parse_quote!{let #ident = &self.#ident;});
			}
			for ident in change_cond_idents.iter().filter(has_ident) {
				change_block.stmts.insert(0, syn::parse_quote!{let #ident = self.#ident;});
			}
			for ident in value_idents.iter().filter(has_ident) {
				value_block.stmts.insert(0, syn::parse_quote!{let #ident = self.#ident;});
			}
			for ident in value_cond_idents.iter().filter(has_ident) {
				value_block.stmts.insert(0, syn::parse_quote!{let #ident = self.#ident;});
			}
			
			(
				kind_type,
				value_block,
				value_idents,
				change_block,
				change_idents,
				flux_type,
				crate_path.to_token_stream(),
			)
		},
		Err(err) => panic!("{}", err),
	};
	
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
		let ident = field.ident.as_ref()
			.expect("identifier must exist");
		
		 // Constant/Linear Types:
		if value_idents.contains(&ident.unraw()) {
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::kind::FluxKind::eval(&#flux::Flux::to_kind(self), time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: moment.#ident,
			};
		}
		
		 // Flux/Moment Types:
		else if change_idents.contains(&ident.unraw()) {
			let field_ty = &field.ty;
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::ToMoment::to_moment(&self.#ident, time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: #flux::Flux::from_moment(moment.#ident),
			};
			field.ty = syn::parse_quote!{
				<#field_ty as #flux::Moment>::Flux
			};
		}
		
		 // Clone Types:
		else {
			moment_fields = quote::quote!{
				#moment_fields
				#ident: self.#ident.clone(),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: moment.#ident,
			};
		}
	}
	
	 // Generate Implementation:
	let ident = item.ident.clone();
	let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();
	let flux_value_impl = quote::quote!{
		#item
		
		impl #impl_generics #flux::Flux for #flux_type #where_clause {
			type Basis = <#kind_type as #flux::Flux>::Basis;
			type Kind = #kind_type;
			
			fn basis(&self) -> Self::Basis
			#value_block
			
			fn change(&self, kind: Self::Kind) -> Self::Kind
			#change_block
		}
		
		impl #impl_generics #flux::ToMoment for #flux_type #where_clause {
			type Moment<'a> = #ident #ty_generics;
			
			fn to_moment(&self, time: #flux::linear::Scalar) -> Self::Moment<'_> {
				#ident { #moment_fields }
			}
		}
		
		impl #impl_generics #flux::ToMomentMut for #flux_type #where_clause {
			type MomentMut<'a> = &'a mut #ident #ty_generics where Self: 'a;
			
			fn to_moment_mut(&mut self, time: #flux::linear::Scalar) -> Self::MomentMut<'_> {
				*self = #flux::ToMoment::to_moment(self, time);
				self
			}
		}
	};
	
	flux_value_impl.into()
}

#[proc_macro_derive(Flux, attributes(flux))]
pub fn derive_flux(item_tokens: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let mut basis_ident: Option<syn::Ident> = None;
	let mut change_closure: Option<syn::ExprClosure> = None;
	
	 // Find `flux` Helper Attribute:
	for field in item.fields.iter() {
		for attr in field.attrs.iter() {
			if attr.meta.path().is_ident("flux") {
				if basis_ident.is_some() {
					panic!("found multiple `flux` helper attributes");
				}
				
				basis_ident = Some(field.ident.clone()
					.expect("basis identifier should exist"));
				
				let syn::Meta::List(meta_list) = &attr.meta
					else { panic!("flux must contain a closure argument") };
				
				change_closure = syn::parse(meta_list.tokens.to_owned().into()).ok();
			}
		}
	}
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
	let basis_ident = basis_ident
		.expect("flux helper attribute should mark a field");
	
	let change_closure = change_closure
		.expect("flux helper attribute should mark a field");
	
	let kind_type = match &change_closure.output {
		syn::ReturnType::Default => syn::parse_quote!{impl #chime::kind::FluxKind},
		syn::ReturnType::Type(_, t) => (**t).clone(),
	};
	
	let flux_impl = quote::quote!{
		impl #impl_params #chime::Flux for #type_name #type_params #impl_clause {
			type Basis = <Self::Kind as #chime::Flux>::Basis;
			type Kind = #kind_type;
			fn basis(&self) -> Self::Basis {
				self.#basis_ident
			}
			fn change(&self, kind_: Self::Kind) -> Self::Kind {
				(#change_closure)(kind_)
			}
		}
	};
	
	flux_impl.into()
}