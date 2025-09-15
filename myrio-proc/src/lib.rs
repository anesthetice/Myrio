use std::ops::Range;

use proc_macro::TokenStream;
use quote::quote;
use syn::{
    Expr, Ident, Type,
    parse::{Parse, ParseStream},
    parse_macro_input,
};

const K_DENSE_VALID_RANGE: Range<usize> = 2..7;

// When `usize` is 64-bit, max{k} = 32 as each nucleotide is repr. by 2 bits
#[cfg(target_pointer_width = "64")]
const K_SPARSE_VALID_RANGE: Range<usize> = 2..33;

// When `usize` is 32-bit, max{k} = 16 as each nucleotide is repr. by 2 bits
#[cfg(target_pointer_width = "32")]
const K_SPARSE_VALID_RANGE: Range<usize> = 2..17;

/// Generates a lengthy match expression required by const generics
/// match k {
///     2 => body!(input, 2),
///     3 => body!(input, 3),
///     ...
///     6 => body!(input, 6),
///     _ => Err(...),
/// }
#[proc_macro]
pub fn gen_match_k_dense(input: TokenStream) -> TokenStream {
    let sequence_expr = parse_macro_input!(input as Expr);

    let arms = K_DENSE_VALID_RANGE.map(|k| {
        quote! {
            #k => body!(#sequence_expr, #k),
        }
    });

    let output = quote! {
        match k {
            #(#arms)*
            _ => panic!("{}", K_DENSE_VALID_RANGE_ERROR_MSG),
        }
    };

    output.into()
}

/// Generates a lengthy match expression required by const generics
/// match k {
///     2 => body!(input, 2),
///     3 => body!(input, 3),
///     ...
///     32 => body!(input, 32),
///     _ => Err(...),
/// }
#[proc_macro]
pub fn gen_match_k_sparse(input: TokenStream) -> TokenStream {
    let sequence_expr = parse_macro_input!(input as Expr);

    let arms = K_SPARSE_VALID_RANGE.map(|k| {
        quote! {
            #k => body!(#sequence_expr, #k),
        }
    });

    let output = quote! {
        match k {
            #(#arms)*
            _ => panic!("{}", MyrSeq::K_SPARSE_VALID_RANGE_ERROR_MSG),
        }
    };

    output.into()
}

struct MacroInput {
    ty: Type,
    op: Ident,
}

impl Parse for MacroInput {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let ty: Type = input.parse()?; // e.g. f64
        let op: Ident = input.parse()?; // e.g. Add
        Ok(MacroInput { ty, op })
    }
}

#[proc_macro]
pub fn impl_ops_for_svec(input: TokenStream) -> TokenStream {
    let MacroInput { ty, op } = parse_macro_input!(input as MacroInput);
    let op_assign = Ident::new(&(op.to_string() + "Assign"), op.span());

    let method = Ident::new(&op.to_string().to_lowercase(), op.span());
    let method_assign = Ident::new(&(method.to_string() + "_assign"), method.span());

    let expanded = quote! {
        impl core::ops::#op<#ty> for SparseVec<#ty> {
            type Output = SparseVec<#ty>;

            fn #method(mut self, rhs: #ty) -> Self::Output {
                self.values_mut().for_each(|v| v.#method_assign(rhs));
                self.sval.#method_assign(rhs);
                self
            }
        }

        impl core::ops::#op<#ty> for &SparseVec<#ty> {
            type Output = SparseVec<#ty>;

            fn #method(self, rhs: #ty) -> Self::Output {
                unsafe {
                    SparseVec::<#ty>::new_unchecked(
                        self.keys.clone(),
                        self.values.iter().map(|v| v.#method(rhs)).collect(),
                        self.dim,
                        self.sval.#method(rhs),
                    )
                }
            }
        }

        impl core::ops::#op_assign<#ty> for SparseVec<#ty> {
            fn #method_assign(&mut self, rhs: #ty) {
                self.values_mut().for_each(|v| v.#method_assign(rhs));
                self.sval.#method_assign(rhs);
            }
        }
    };

    expanded.into()
}
