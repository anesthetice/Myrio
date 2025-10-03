// Imports
use bincode::{BorrowDecode, Decode, Encode};

use crate::tax::clade::Rank;

pub struct Branch<B> {
    pub name: Box<str>,
    pub children: Box<[Node<B>]>,
    pub extra: B,
}

impl<B> Branch<B> {
    pub fn gather_leaves<'a>(&'a self) -> Vec<&'a Leaf> {
        let mut leaves: Vec<&'a Leaf> = Vec::new();
        for node in self.children.iter() {
            node.gather_leaves(&mut leaves);
        }
        leaves
    }
}

#[derive(Debug, Clone, Encode, Decode)]
pub struct Leaf {
    pub name: Box<str>,
    pub payload_id: usize,
}

pub enum Node<B> {
    Branch(Branch<B>),
    Leaf(Leaf),
}

impl<B> Node<B> {
    pub fn new_branch(
        name: Box<str>,
        children: Box<[Node<B>]>,
        extra: B,
    ) -> Self {
        Self::Branch(Branch { name, children, extra })
    }

    pub fn new_leaf(
        name: Box<str>,
        payload_id: usize,
    ) -> Self {
        Self::Leaf(Leaf { name, payload_id })
    }

    pub fn get_name(&self) -> &str {
        match self {
            Self::Branch(branch) => &branch.name,
            Self::Leaf(leaf) => &leaf.name,
        }
    }

    pub fn gather_leaves<'a>(
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

    pub fn gather_leaves_mut<'a>(
        &'a mut self,
        output: &mut Vec<&'a mut Leaf>,
    ) {
        match self {
            Self::Branch(branch) => {
                for node in branch.children.iter_mut() {
                    node.gather_leaves_mut(output);
                }
            }
            Self::Leaf(leaf) => output.push(leaf),
        }
    }

    pub fn gather_branches<'a>(
        &'a self,
        output: &mut Vec<&'a Branch<B>>,
    ) {
        if let Self::Branch(branch) = self {
            output.push(branch);
            branch.children.iter().for_each(|node| node.gather_branches(output));
        }
    }

    pub fn gather_leaves_and_branches<'a>(
        &'a self,
        output_leaves: &mut Vec<&'a Leaf>,
        output_branches: &mut Vec<&'a Branch<B>>,
    ) {
        match self {
            Self::Branch(branch) => {
                output_branches.push(branch);
                branch
                    .children
                    .iter()
                    .for_each(|node| node.gather_leaves_and_branches(output_leaves, output_branches));
            }
            Self::Leaf(leaf) => output_leaves.push(leaf),
        }
    }

    pub fn unwrap_branch(self) -> Branch<B> {
        match self {
            Self::Branch(branch) => branch,
            _ => unreachable!(),
        }
    }

    pub fn unwrap_branch_ref(&self) -> &Branch<B> {
        match self {
            Self::Branch(branch) => branch,
            _ => unreachable!(),
        }
    }

    pub fn unwrap_leaf(self) -> Leaf {
        match self {
            Self::Leaf(leaf) => leaf,
            _ => unreachable!(),
        }
    }

    pub fn unwrap_leaf_ref(&self) -> &Leaf {
        match self {
            Self::Leaf(leaf) => leaf,
            _ => unreachable!(),
        }
    }
}

impl<B: Clone> std::clone::Clone for Node<B> {
    fn clone(&self) -> Self {
        match &self {
            Self::Branch(branch) => {
                Self::new_branch(branch.name.clone(), branch.children.clone(), branch.extra.clone())
            }
            Self::Leaf(leaf) => Self::Leaf(leaf.clone()),
        }
    }
}

pub struct TaxTreeCore<B, L> {
    pub gene: String,
    pub root: Node<B>,
    pub root_rank: Rank,
    pub payloads: Box<[L]>,
}

impl<B, L> TaxTreeCore<B, L> {
    pub fn new(
        gene: String,
        root: Node<B>,
        root_rank: Rank,
        payloads: Box<[L]>,
    ) -> Self {
        Self { gene, root, root_rank, payloads }
    }

    pub fn gather_leaves(&self) -> Vec<&Leaf> {
        let mut output = Vec::with_capacity(self.payloads.len());
        self.root.gather_leaves(&mut output);
        output
    }

    pub fn gather_leaves_mut(&mut self) -> Vec<&mut Leaf> {
        let mut output = Vec::with_capacity(self.payloads.len());
        self.root.gather_leaves_mut(&mut output);
        output
    }

    pub fn gather_branches_at_rank(
        &self,
        rank: Rank,
    ) -> Vec<&Branch<B>> {
        if rank > self.root_rank {
            return Vec::with_capacity(0);
        }

        fn dive_recursive<'a, B>(
            node: &'a Node<B>,
            current_height: usize,
            desired_height: usize,
            output: &mut Vec<&'a Branch<B>>,
        ) {
            if let Node::Branch(branch) = node {
                if current_height == desired_height {
                    output.push(branch);
                } else {
                    for child in branch.children.iter() {
                        dive_recursive(child, current_height - 1, desired_height, output);
                    }
                }
            }
        }

        let mut output = Vec::new();
        dive_recursive(&self.root, self.root_rank as usize, rank as usize, &mut output);
        output
    }

    pub fn gather_branches_mut_at_rank(
        &mut self,
        rank: Rank,
    ) -> Vec<&mut Branch<B>> {
        if rank > self.root_rank {
            return Vec::with_capacity(0);
        }

        fn dive_recursive<'a, B>(
            node: &'a mut Node<B>,
            current_height: usize,
            desired_height: usize,
            output: &mut Vec<&'a mut Branch<B>>,
        ) {
            if let Node::Branch(branch) = node {
                if current_height == desired_height {
                    output.push(branch);
                } else {
                    for child in branch.children.iter_mut() {
                        dive_recursive(child, current_height - 1, desired_height, output);
                    }
                }
            }
        }

        let mut output = Vec::new();
        dive_recursive(&mut self.root, self.root_rank as usize, rank as usize, &mut output);
        output
    }

    pub fn excise(
        &mut self,
        name: &str,
        rank: Rank,
    ) {
        // rename to target_name for clarity
        let target_name = name;

        let mut payload_ids_to_excise: Vec<usize> = Vec::new();

        let branches_above = self.gather_branches_mut_at_rank(rank.above().unwrap());
        for branch in branches_above {
            if branch.children.iter().all(|node| node.get_name() != target_name) {
                continue;
            }
            let mut children = std::mem::take(&mut branch.children).into_vec();

            children.retain(|node| {
                if node.get_name() != target_name {
                    true
                } else {
                    let mut leaves = Vec::new();
                    node.gather_leaves(&mut leaves);
                    payload_ids_to_excise.extend(leaves.iter().map(|leaf| leaf.payload_id));
                    false
                }
            });

            branch.children = children.into_boxed_slice();
        }

        self.remove_payloads(&payload_ids_to_excise);

        #[cfg(debug_assertions)]
        println!("Excision for the {} tree affected {} leaves", self.gene, payload_ids_to_excise.len());
    }

    /// Method used to remove multiple payloads by their associated id, assumes associated leaves/branches have already been removed, will panic if such leaves are encountered
    fn remove_payloads(
        &mut self,
        payload_ids: &[usize],
    ) {
        if payload_ids.is_empty() {
            return;
        }

        // Extract payloads from self
        let mut payloads = std::mem::take(&mut self.payloads).into_vec();

        // Sort and dedup `payload_ids` (important for step 2)
        let mut payload_ids = payload_ids.to_vec();
        payload_ids.sort_unstable();
        payload_ids.dedup();

        // Step 1, remove all payloads whose index matches one specified in `payload_ids`
        {
            let mut should_be_retained_vec: Vec<bool> = vec![true; payloads.len()];
            for pid in payload_ids.iter().copied() {
                should_be_retained_vec[pid] = false;
            }
            let mut should_be_retained_iter = should_be_retained_vec.into_iter();
            payloads.retain(|_| should_be_retained_iter.next().unwrap());
        }

        // Step 2, adjust the `payload_id` of every leaf
        for leaf in self.gather_leaves_mut() {
            let Err(left_shift) = payload_ids.binary_search(&leaf.payload_id) else {
                std::panic!("Found an unexpected leaf: {leaf:?}");
            };
            leaf.payload_id -= left_shift;
        }

        // Inject payloads back into self
        self.payloads = payloads.into_boxed_slice()
    }
}

impl<B, L> core::clone::Clone for TaxTreeCore<B, L>
where
    B: Clone,
    L: Clone,
{
    fn clone(&self) -> Self {
        Self {
            gene: self.gene.clone(),
            root: self.root.clone(),
            root_rank: self.root_rank,
            payloads: self.payloads.clone(),
        }
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
        bincode::Encode::encode(&self.root, encoder)?;
        bincode::Encode::encode(&self.root_rank, encoder)?;
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
            root: bincode::Decode::decode(decoder)?,
            root_rank: bincode::Decode::decode(decoder)?,
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
            root: bincode::BorrowDecode::borrow_decode(decoder)?,
            root_rank: bincode::BorrowDecode::borrow_decode(decoder)?,
            payloads: bincode::BorrowDecode::borrow_decode(decoder)?,
        })
    }
}
