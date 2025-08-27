// Imports
use crate::tax::clade::Rank;
use bincode::{BorrowDecode, Decode, Encode};

pub(super) struct Branch<B> {
    pub(super) name: Box<str>,
    pub(super) children: Box<[Node<B>]>,
    pub(super) extra: B,
}

#[derive(Encode, Decode)]
pub(super) struct Leaf {
    pub(super) name: Box<str>,
    pub(super) payload_id: usize,
}

pub(super) enum Node<B> {
    Branch(Branch<B>),
    Leaf(Leaf),
}

pub(super) struct TaxTreeCore<B, L> {
    pub(super) gene: String,
    pub(super) highest_rank: Rank,
    pub(super) roots: Box<[Node<B>]>,
    pub(super) payloads: Box<[L]>,
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
