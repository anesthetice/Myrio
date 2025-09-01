// Imports
use bincode::{BorrowDecode, Decode, Encode};

use crate::tax::clade::Rank;

pub(crate) struct Branch<B> {
    pub(crate) name: Box<str>,
    pub(crate) children: Box<[Node<B>]>,
    pub(crate) extra: B,
}

#[derive(Clone, Encode, Decode)]
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

    pub(crate) fn gather_leaves<'a>(
        &'a self,
        output: &mut Vec<&'a Leaf>,
    ) {
        match self {
            Self::Branch(branch) => {
                for node in &branch.children {
                    node.gather_leaves(output);
                }
            }
            Self::Leaf(leaf) => output.push(leaf),
        }
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
        let mut output = Vec::new();
        for root in self.roots.iter() {
            root.gather_leaves(&mut output);
        }
        output
    }
}

impl<B, L> std::fmt::Display for TaxTreeCore<B, L>
where
    B: std::fmt::Display,
    L: std::fmt::Display,
{
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        writeln!(f, "Root ({})", self.gene)?;

        fn space(
            depth: u8,
            is_end_stack: u8,
        ) -> String {
            let mut s = String::new();
            for i in 0..depth {
                if is_end_stack & (0b1_u8 << i) == (0b1_u8 << i) {
                    s.push_str("    ");
                } else {
                    s.push_str("│   ");
                }
            }
            s
        }

        fn dive<B, L>(
            node: &Node<B>,
            payloads: &[L],
            depth: u8,
            is_end: bool,
            mut is_end_stack: u8,
            f: &mut std::fmt::Formatter<'_>,
        ) -> std::result::Result<(), std::fmt::Error>
        where
            B: std::fmt::Display,
            L: std::fmt::Display,
        {
            match node {
                Node::Branch(branch) => {
                    writeln!(
                        f,
                        "{}{}── {} ({})",
                        space(depth, is_end_stack),
                        if is_end { "└" } else { "├" },
                        branch.name,
                        branch.extra
                    )?;
                    if is_end {
                        is_end_stack |= 0b1_u8 << depth
                    }
                    let len = branch.children.len();
                    match len {
                        0 => (),
                        1 => dive(&branch.children[0], payloads, depth + 1, true, is_end_stack, f)?,
                        2.. => {
                            for node in branch.children[0..len - 1].iter() {
                                dive(node, payloads, depth + 1, false, is_end_stack, f)?;
                            }
                            dive(
                                branch.children.last().unwrap(),
                                payloads,
                                depth + 1,
                                true,
                                is_end_stack,
                                f,
                            )?;
                        }
                    }
                }
                Node::Leaf(leaf) => {
                    writeln!(
                        f,
                        "{}{}── {} ({})",
                        space(depth, is_end_stack),
                        if is_end { "└" } else { "├" },
                        leaf.name,
                        payloads[leaf.payload_id],
                    )?;
                }
            }
            Ok(())
        }

        let len = self.roots.len();
        match len {
            0 => (),
            1 => dive(&self.roots[0], &self.payloads, 0, true, 0b0, f)?,
            2.. => {
                for node in self.roots[0..len - 1].iter() {
                    dive(node, &self.payloads, 0, false, 0b0, f)?;
                }
                dive(self.roots.last().unwrap(), &self.payloads, 0, true, 0b0, f)?;
            }
        }

        Ok(())
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

#[cfg(test)]
mod test {
    use indoc::indoc;

    use super::*;

    #[test]
    fn display_test() {
        #[rustfmt::skip]
        let tree: TaxTreeCore<f64, f64> = TaxTreeCore {
            gene: "Test".to_string(),
            highest_rank: Rank::Family,
            roots: [
                Node::new_branch(
                    "test-1".into(),
                    [
                        Node::new_branch(
                            "test-11".into(),
                            [
                                Node::new_leaf("test-111".into(), 0),
                                Node::new_leaf("test-112".into(), 1),
                            ].into(),
                            0.15,
                        ),
                        Node::new_branch(
                            "test-12".into(),
                            [
                                Node::new_leaf("test-121".into(), 2),
                                Node::new_leaf("test-122".into(), 3),
                            ].into(),
                            0.0,
                        )
                        ].into(),
                    0.075),
                Node::new_branch(
                    "test-2".into(),
                    [
                        Node::new_branch(
                            "test-21".into(),
                            [
                                Node::new_leaf("test-211".into(), 4),
                            ].into(),
                            0.9,
                        )
                        ].into(),
                    0.9),
                ].into(),
            payloads: [0.1, 0.2, 0.0, 0.0, 0.9].into()
        };

        eprintln!("{}", tree);

        let expected = indoc! {"
            Root (Test)
            ├── test-1 (0.075)
            │   ├── test-11 (0.15)
            │   │   ├── test-111 (0.1)
            │   │   └── test-112 (0.2)
            │   └── test-12 (0)
            │       ├── test-121 (0)
            │       └── test-122 (0)
            └── test-2 (0.9)
                └── test-21 (0.9)
                    └── test-211 (0.9)
        "};
        assert_eq!(tree.to_string().as_str(), expected);
    }
}
