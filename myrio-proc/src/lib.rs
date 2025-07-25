use std::ops::Range;

use proc_macro::TokenStream;
use quote::quote;
use syn::{Expr, parse_macro_input};

const K_VALID_RANGE: Range<usize> = 2..43;

/// Generates a lengthy match expression required by const generics
/// match k {
///     2 => body!(self.sequence, 2),
///     3 => body!(self.sequence, 3),
///     ...
///     42 => body!(self.sequence, 42),
///     _ => Err(Error::InvalidKmerSize),
/// }
#[proc_macro]
pub fn match_k(input: TokenStream) -> TokenStream {
    let sequence_expr = parse_macro_input!(input as Expr);

    let arms = K_VALID_RANGE.map(|k| {
        quote! {
            #k => body!(#sequence_expr, #k),
        }
    });

    let output = quote! {
        match k {
            #(#arms)*
            _ => Err(Error::InvalidKmerSize),
        }
    };

    output.into()
}
