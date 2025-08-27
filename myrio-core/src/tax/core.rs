// Imports
use crate::tax::clade::Rank;
use bincode::{BorrowDecode, Decode, Encode};

pub(crate) struct Branch<B> {
    pub(crate) name: Box<str>,
    pub(crate) children: Box<[Node<B>]>,
    pub(crate) extra: B,
}

#[derive(Encode, Decode)]
pub(crate) struct Leaf {
    pub(crate) name: Box<str>,
    pub(crate) payload_id: usize,
}

pub(crate) enum Node<B> {
    Branch(Branch<B>),
    Leaf(Leaf),
}

impl<B> Node<B> {
    pub(crate) fn new_branch(
        name: Box<str>,
        children: Box<[Node<B>]>,
        extra: B,
    ) -> Self {
        Self::Branch(Branch { name, children, extra })
    }

    pub(crate) fn new_leaf(
        name: Box<str>,
        payload_id: usize,
    ) -> Self {
        Self::Leaf(Leaf { name, payload_id })
    }
}

pub(crate) struct TaxTreeCore<B, L> {
    pub(crate) gene: String,
    pub(crate) highest_rank: Rank,
    pub(crate) roots: Box<[Node<B>]>,
    pub(crate) payloads: Box<[L]>,
}

impl<B, L> TaxTreeCore<B, L> {
    pub(crate) fn new(
        gene: String,
        highest_rank: Rank,
        roots: Box<[Node<B>]>,
        payloads: Box<[L]>,
    ) -> Self {
        Self { gene, highest_rank, roots, payloads }
    }
    pub fn gather_leaves(&self) -> Vec<&Leaf> {
        fn recursive_dive<'a, B>(
            node: &'a Node<B>,
            output: &mut Vec<&'a Leaf>,
        ) {
            match node {
                Node::Branch(branch) => {
                    for node in &branch.children {
                        recursive_dive(node, output);
                    }
                }
                Node::Leaf(leaf) => output.push(leaf),
            }
        }
        let mut output = Vec::new();
        for root in self.roots.iter() {
            recursive_dive(root, &mut output);
        }
        output
    }
}

impl<B: Encode> bincode::Encode for Branch<B> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.name, encoder)?;
        bincode::Encode::encode(&self.children, encoder)?;
        bincode::Encode::encode(&self.extra, encoder)?;
        Ok(())
    }
}

impl<Context, B: Decode<Context> + 'static> bincode::Decode<Context> for Branch<B> {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            name: bincode::Decode::decode(decoder)?,
            children: bincode::Decode::decode(decoder)?,
            extra: bincode::Decode::decode(decoder)?,
        })
    }
}

impl<'de, Context, B: BorrowDecode<'de, Context> + 'de> bincode::BorrowDecode<'de, Context> for Branch<B> {
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            name: bincode::BorrowDecode::borrow_decode(decoder)?,
            children: bincode::BorrowDecode::borrow_decode(decoder)?,
            extra: bincode::BorrowDecode::borrow_decode(decoder)?,
        })
    }
}

impl<B: Encode> bincode::Encode for Node<B> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        match self {
            Node::Branch(b) => {
                0u8.encode(encoder)?;
                b.encode(encoder)
            }
            Node::Leaf(l) => {
                1u8.encode(encoder)?;
                l.encode(encoder)
            }
        }
    }
}

impl<Context, B: Decode<Context> + 'static> bincode::Decode<Context> for Node<B> {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        let tag: u8 = Decode::decode(decoder)?;
        match tag {
            0 => Ok(Node::Branch(Decode::decode(decoder)?)),
            1 => Ok(Node::Leaf(Decode::decode(decoder)?)),
            _ => Err(bincode::error::DecodeError::OtherString("invalid enum tag".into())),
        }
    }
}

impl<'de, Context, B: BorrowDecode<'de, Context> + 'de> bincode::BorrowDecode<'de, Context> for Node<B> {
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        let tag: u8 = BorrowDecode::borrow_decode(decoder)?;
        match tag {
            0 => Ok(Node::Branch(BorrowDecode::borrow_decode(decoder)?)),
            1 => Ok(Node::Leaf(BorrowDecode::borrow_decode(decoder)?)),
            _ => Err(bincode::error::DecodeError::OtherString("invalid enum tag".into())),
        }
    }
}

impl<B: Encode, L: Encode> bincode::Encode for TaxTreeCore<B, L> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.gene, encoder)?;
        bincode::Encode::encode(&self.highest_rank, encoder)?;
        bincode::Encode::encode(&self.roots, encoder)?;
        bincode::Encode::encode(&self.payloads, encoder)?;
        Ok(())
    }
}

impl<Context, B, L> bincode::Decode<Context> for TaxTreeCore<B, L>
where
    B: Decode<Context> + 'static,
    L: Decode<Context> + 'static,
{
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            gene: bincode::Decode::decode(decoder)?,
            highest_rank: bincode::Decode::decode(decoder)?,
            roots: bincode::Decode::decode(decoder)?,
            payloads: bincode::Decode::decode(decoder)?,
        })
    }
}

impl<'de, Context, B, L> bincode::BorrowDecode<'de, Context> for TaxTreeCore<B, L>
where
    B: BorrowDecode<'de, Context> + 'de,
    L: BorrowDecode<'de, Context> + 'de,
{
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            gene: bincode::BorrowDecode::borrow_decode(decoder)?,
            highest_rank: bincode::BorrowDecode::borrow_decode(decoder)?,
            roots: bincode::BorrowDecode::borrow_decode(decoder)?,
            payloads: bincode::BorrowDecode::borrow_decode(decoder)?,
        })
    }
}
