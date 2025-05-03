use proc_macro::*;
use quote::ToTokens;

static BASIS_HELP: &'static str = "\
	\nhelp: mark a single field as the basis by placing `#[basis]` above it\
";

static CHANGE_HELP: &'static str = "\
	\nhelp: `#[change(op)]` is used to mark fields that apply an operation over time:\
	\n        `add_per(u)` - Add per unit of `chime::time::Time` (fields of type `impl Flux`)\
	\n        `sub_per(u)` - Subtract per unit of `chime::time::Time` (fields of type `impl Flux`)\
	\n        `add`        - Add directly (fields of type `impl Poly` or `Change<impl Flux>`)\
	\n        `sub`        - Subtract directly (fields of type `impl Poly` or `Change<impl Flux>`)\
";

static MOMENT_HELP: &'static str = "";

#[proc_macro_derive(Flux, attributes(basis, change))]
pub fn derive_flux(item_tokens: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let mut basis: Option<(syn::Member, syn::Type)> = None;
	let mut change_expr: syn::Expr = syn::parse_quote!{#chime::Flux::accum(self).into_change()};
	let mut change_type: syn::Type = syn::parse_quote!{#chime::change::Nil<Self::Basis>};
	
	 // Find Helper Attributes:
	for (field_index, field) in item.fields.iter().enumerate() {
		let field_member = field.ident.clone()
			.map(Into::<syn::Member>::into)
			.unwrap_or_else(|| field_index.into());
		
		let field_type = &field.ty;
		
		for attr in field.attrs.iter() {
			if attr.meta.path().is_ident("basis") {
				if basis.is_some() {
					panic!("multiple fields declared as basis{}", BASIS_HELP);
				}
				basis = Some((field_member.clone(), field_type.clone()));
				continue
			}
			if attr.meta.path().is_ident("change") {
				let syn::Meta::List(meta_list) = &attr.meta
					else { panic!("expected argument list for `{}`{}",
						attr.to_token_stream(), CHANGE_HELP) };
				
				match meta_list.parse_args::<syn::Expr>() {
					Ok(syn::Expr::Path(syn::ExprPath { path, .. })) => {
						let mut rhs: syn::Expr = syn::parse_quote!{self.#field_member};
						
						let _op: syn::LitChar = loop {
							 // Apply by Copy:
							if path.is_ident("add") { break syn::parse_quote!{'+'} }
							if path.is_ident("sub") { break syn::parse_quote!{'-'} }
							if path.is_ident("mul") { break syn::parse_quote!{'*'} }
							if path.is_ident("div") { break syn::parse_quote!{'/'} }

							 // Apply by Reference:
							rhs = syn::parse_quote!{&#rhs};
							if path.is_ident("add_ref") { break syn::parse_quote!{'+'} }
							if path.is_ident("sub_ref") { break syn::parse_quote!{'-'} }
							if path.is_ident("mul_ref") { break syn::parse_quote!{'*'} }
							if path.is_ident("div_ref") { break syn::parse_quote!{'/'} }

							panic!("invalid change operation, `{}`{}",
								path.to_token_stream(), CHANGE_HELP)
						};
						
						change_expr = syn::parse_quote!{
							(#change_expr >> #rhs)
						};
						change_type = syn::parse_quote!{
							<#change_type as std::ops::Shr<#field_type>>::Output
						};
					},
					Ok(syn::Expr::Call(syn::ExprCall { func, args, .. })) => {
						let syn::Expr::Path(syn::ExprPath { path, .. }) = &*func
							else { panic!("invalid change operation, `{}`{}",
								func.to_token_stream(), CHANGE_HELP) };
						
						let (op, op_trait): (syn::BinOp, syn::TypePath) = loop {
							if path.is_ident("add_per") { break (syn::parse_quote!{+}, syn::parse_quote!{std::ops::Add}) }
							if path.is_ident("sub_per") { break (syn::parse_quote!{-}, syn::parse_quote!{std::ops::Sub}) }
							if path.is_ident("mul_per") { break (syn::parse_quote!{*}, syn::parse_quote!{std::ops::Mul}) }
							if path.is_ident("div_per") { break (syn::parse_quote!{/}, syn::parse_quote!{std::ops::Div}) }
							
							panic!("invalid change operation, `{}`{}",
								path.to_token_stream(), CHANGE_HELP)
						};
						
						 // Get Unit:
						if args.len() != 1 {
							panic!("expected one argument for `{}`{}",
								meta_list.tokens, CHANGE_HELP);
						}
						let unit = args.first().unwrap();
						
						change_expr = syn::parse_quote!{
							(#change_expr #op #chime::Flux::per(&self.#field_member, #unit))
						};
						change_type = syn::parse_quote!{
							<#change_type as #op_trait<#chime::Rate<#field_type>>>::Output
						};
					},
					Ok(meta) => panic!("invalid change operation, `{}`{}",
						meta.to_token_stream(), CHANGE_HELP),
					Err(e) => panic!("{}{}", e, CHANGE_HELP),
				}
				continue
			}
		}
	}
	
	let (basis_member, basis_type) = basis.unwrap_or_else(|| panic!(
		"no basis declared{}",
		BASIS_HELP
	));
	
	let trait_impl = quote::quote!{
		impl #impl_params #chime::Flux for #type_name #type_params #impl_clause {
			type Basis = #basis_type;
			type Change = #change_type;
			fn basis(&self) -> Self::Basis {
				self.#basis_member
			}
			fn change(&self) -> Self::Change {
				#change_expr
			}
		}
	};
	
	trait_impl.into()
}

#[proc_macro_derive(ToMoment, attributes(moment, basis))]
pub fn derive_to_moment(item_tokens: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let mut fields = syn::punctuated::Punctuated::new();
	
	 // Build Struct Literal:
	for (field_index, field) in item.fields.iter().enumerate() {
		let member = field.ident.clone()
			.map(Into::<syn::Member>::into)
			.unwrap_or_else(|| field_index.into());
		
		let mut expr = syn::parse_quote!{
			#chime::ToMoment::to_moment(&self.#member, time)
		};
		
		 // Helper Attributes:
		for attr in field.attrs.iter() {
			if attr.meta.path().is_ident("basis") {
				expr = syn::parse_quote!{
					#chime::poly::Poly::eval(
						&#chime::Flux::to_poly(self),
						#chime::linear::Linear::from_f64(time),
					)
				};
				continue
			}
			if attr.meta.path().is_ident("moment") {
				let syn::Meta::List(meta_list) = &attr.meta
					else { panic!("expected argument list for `{}`{}",
						attr.to_token_stream(), MOMENT_HELP) };
				
				let Ok(arg_ident) = meta_list.parse_args::<syn::Ident>()
					else { panic!("expected `moment(ignore)`, got `{}`{}",
						attr.to_token_stream(), MOMENT_HELP) };
				
				if arg_ident != "ignore" {
					panic!("expected `moment(ignore)`, got `{}`{}",
						attr.to_token_stream(), MOMENT_HELP)
				}
				
				expr = syn::parse_quote!{self.#member.clone()};
				continue
			}
		}
		
		fields.push(syn::FieldValue {
			attrs: vec![],
			member,
			colon_token: Some(Default::default()),
			expr,
		});
	}
	
	let moment_struct = syn::ExprStruct {
        attrs: vec![],
        qself: None,
        path: syn::parse_quote!{Self},
        brace_token: Default::default(),
        fields,
        dot2_token: None,
        rest: None,
	};
	
	let trait_impl = quote::quote!{
		impl #impl_params #chime::ToMoment for #type_name #type_params #impl_clause {
			type Moment<'a_> = Self; // !!! Replace marked type params later
			fn to_moment(&self, time: f64) -> Self::Moment<'_> {
				#moment_struct
			}
		}
	};
	
	trait_impl.into()
}

#[proc_macro_derive(ToMomentMut)]
pub fn derive_to_moment_mut(item_tokens: TokenStream) -> TokenStream {
	let item: syn::ItemStruct = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let trait_impl = quote::quote!{
		impl #impl_params #chime::ToMomentMut for #type_name #type_params #impl_clause {
			type MomentMut<'a_> = &'a_ mut Self;
			fn to_moment_mut(&mut self, time: f64) -> Self::MomentMut<'_> {
				*self = #chime::ToMoment::to_moment(&self, time);
				self
			}
		}
	};
	
	trait_impl.into()
}

#[proc_macro_derive(TemporalComponent)]
pub fn temporal_component_derive(item_tokens: TokenStream) -> TokenStream {
	let mut item: syn::DeriveInput = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
    item.generics
        .make_where_clause()
        .predicates
        .push(syn::parse_quote!{Self: Send + Sync + 'static});
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let trait_impl = quote::quote!{
		impl #impl_params #chime::temporal::TemporalComponent for #type_name #type_params #impl_clause {
		}
	};
	
	trait_impl.into()
}

#[proc_macro_derive(TemporalResource)]
pub fn temporal_resource_derive(item_tokens: TokenStream) -> TokenStream {
	let mut item: syn::DeriveInput = match syn::parse(item_tokens) {
		Ok(item) => item,
		Err(e) => panic!("{}", e),
	};
	
	let chime: syn::Path = syn::parse_quote!{chime};
	
    item.generics
        .make_where_clause()
        .predicates
        .push(syn::parse_quote!{Self: Send + Sync + 'static});
	
	let type_name = item.ident.clone();
	let (impl_params, type_params, impl_clause) = item.generics.split_for_impl();
	
	let trait_impl = quote::quote!{
		impl #impl_params #chime::temporal::TemporalResource for #type_name #type_params #impl_clause {
		}
	};
	
	trait_impl.into()
}