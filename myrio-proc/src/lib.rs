use std::ops::Range;

use proc_macro::TokenStream;
use quote::quote;
use syn::{Expr, Ident, parse_macro_input};

const K_DENSE_VALID_RANGE: Range<usize> = 2..10;
const K_SPARSE_VALID_RANGE: Range<usize> = 2..43;

/// Generates a lengthy match expression required by const generics
/// match k {
///     2 => body!(self.sequence, 2),
///     3 => body!(self.sequence, 3),
///     ...
///     9 => body!(self.sequence, 9),
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
            _ => Err(Error::InvalidKmerSize(Self::K_DENSE_VALID_RANGE_ERROR_MSG)),
        }
    };

    output.into()
}

/// Generates a lengthy match expression required by const generics
/// match k {
///     2 => body!(self.sequence, 2),
///     3 => body!(self.sequence, 3),
///     ...
///     42 => body!(self.sequence, 42),
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
            _ => Err(Error::InvalidKmerSize(Self::K_SPARSE_VALID_RANGE_ERROR_MSG)),
        }
    };

    output.into()
}

#[proc_macro]
pub fn impl_f64_ops_for_sfvec(input: TokenStream) -> TokenStream {
    let op = parse_macro_input!(input as Ident);
    let op_assign = Ident::new(&(op.to_string() + "Assign"), op.span());

    let method = Ident::new(&op.to_string().to_lowercase(), op.span());
    let method_assign = Ident::new(&(method.to_string() + "_assign"), method.span());

    let expanded = quote! {
        impl core::ops::#op<f64> for SparseFloatVec {
            type Output = SparseFloatVec;

            fn #method(mut self, rhs: f64) -> Self::Output {
                self.values_mut().for_each(|v| v.#method_assign(rhs));
                self.sval.#method_assign(rhs);
                self
            }
        }

        impl core::ops::#op<f64> for &SparseFloatVec {
            type Output = SparseFloatVec;

            fn #method(self, rhs: f64) -> Self::Output {
                unsafe {
                    SparseFloatVec::new_unchecked(
                        self.keys.clone(),
                        self.values.iter().map(|v| v.#method(rhs)).collect(),
                        self.dim,
                        self.sval.#method(rhs),
                    )
                }
            }
        }

        impl core::ops::#op_assign<f64> for SparseFloatVec {
            fn #method_assign(&mut self, rhs: f64) {
                self.values_mut().for_each(|v| v.#method_assign(rhs));
                self.sval.#method_assign(rhs);
            }
        }
    };

    expanded.into()
}
