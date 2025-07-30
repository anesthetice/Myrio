use std::ops::Range;

use proc_macro::TokenStream;
use quote::quote;
use syn::{Expr, parse_macro_input};

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
