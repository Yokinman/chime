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
			
			let mut value_block: syn::Block = syn::parse_quote!{{ #value_expr }};
			let mut change_block: syn::Block = syn::parse_quote!{{ (#change_expr)(accum) }};
			
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
	let mut moment_value_type: syn::Type = syn::parse_quote!{
		<<Self::Flux as #flux::Flux>::Kind as #flux::kind::FluxKind>::Value
	};
	let mut moment_fields = Default::default();
	let mut flux_fields = Default::default();
	for field in &mut flux_item.fields {
		let ident = field.ident.as_ref()
			.expect("identifier must exist");
		
		 // Constant/Linear Types:
		if value_idents.contains(&ident.unraw()) {
			let field_ty = syn::parse_quote!{
				<<Self as #flux::_hidden::InnerFlux>::Kind
				as #flux::kind::FluxKind>::Value
			};
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::Flux::value(&#flux::FluxRef::new(&self, base_time), time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: self.#ident,
			};
			moment_value_type = std::mem::replace(&mut field.ty, field_ty);
		}
		
		 // Flux/Moment Types:
		else if change_idents.contains(&ident.unraw()) {
			let field_ty = &field.ty;
			field.ty = syn::parse_quote!{
				<#field_ty as #flux::Moment>::Flux
			};
			moment_fields = quote::quote!{
				#moment_fields
				#ident: #flux::Flux::to_moment(#flux::FluxValue::new(self.#ident, base_time), time),
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: #flux::Flux::into_inner_flux(#flux::Moment::to_flux(self.#ident, time)),
			};
		}
		
		 // Copy Types:
		else {
			moment_fields = quote::quote!{
				#moment_fields
				#ident: self.#ident,
			};
			flux_fields = quote::quote!{
				#flux_fields
				#ident: self.#ident,
			};
		}
	}
	
	 // Generate Implementation:
	let ident = item.ident.clone();
	let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();
	let flux_value_impl = quote::quote!{
		#item
		
		impl #impl_generics #flux::Moment for #ident #ty_generics #where_clause {
			type Flux = #flux::FluxValue<#flux_type>;
			type Value = <#moment_value_type as #flux::linear::LinearPlus>::Inner;
			fn to_flux(self, time: #flux::time::Time) -> Self::Flux {
				#flux::FluxValue::new(#flux_type { #flux_fields }, time)
			}
		}
		
		impl #impl_generics #flux::_hidden::InnerFlux for #flux_type #where_clause {
			type Moment = #ident #ty_generics;
			type Kind = #kind_type;
			
			fn base_value(&self, base_time: #flux::time::Time) -> <Self::Kind as #flux::kind::FluxKind>::Value
			#value_block
			
			fn change<'a>(&self, accum: <Self::Kind as #flux::kind::FluxKind>::Accum<'a>)
				-> <Self::Kind as #flux::kind::FluxKind>::OutAccum<'a>
			#change_block
			
			fn to_moment(self, base_time: #flux::time::Time, time: #flux::time::Time) -> Self::Moment {
				#ident { #moment_fields }
			}
		}
	};
	
	flux_value_impl.into()
}