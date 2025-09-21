use console::style;

use crate::tax::{
    clade::Rank,
    core::{Node, TaxTreeCore},
};

impl<B, L> std::fmt::Display for TaxTreeCore<B, L>
where
    B: std::fmt::Display,
    L: std::fmt::Display,
{
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        let root_branch = self.root.unwrap_branch_ref();
        let rbc_len = root_branch.children.len();

        let title = format!("-- {} taxonomic tree --", self.gene);
        writeln!(f, "{}", style(title).bold())?;

        if rbc_len == 0 {
            return writeln!(f, "()");
        }

        fn to_dim_with_parenthesis(rank: Rank) -> console::StyledObject<String> {
            style(format!("({})", rank.to_char())).dim()
        }

        writeln!(
            f,
            "{} {} {}",
            root_branch.name,
            to_dim_with_parenthesis(self.root_rank),
            root_branch.extra,
        )?;

        #[rustfmt::skip]
        fn dive_recursive<B, L>(
            node: &Node<B>,
            payloads: &[L],
            is_end: bool,
            mut left: String,
            rank: Rank,
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
                        "{}┃\n{}{}━━━━━━ {} {} {}",
                        left,
                        left,
                        if is_end { "┗" } else { "┣" },
                        branch.name,
                        to_dim_with_parenthesis(rank),
                        branch.extra,
                    )?;

                    left += if is_end {"        "} else {"┃       "};

                    if !branch.children.is_empty() {
                        let len = branch.children.len();
                        for node in branch.children[0..len - 1].iter() {
                            dive_recursive(node, payloads, false, left.clone(), rank.below().unwrap(), f)?;
                        }
                        dive_recursive(&branch.children[len - 1], payloads, true, left.clone(), rank.below().unwrap(), f)?;
                    }
                }
                Node::Leaf(leaf) => {
                    writeln!(
                        f,
                        "{}{}─── {} {} {}",
                        left,
                        if is_end { "└" } else { "├" },
                        leaf.name,
                        to_dim_with_parenthesis(rank),
                        payloads[leaf.payload_id],
                    )?;
                }
            }
            Ok(())
        }
        let rank_below = self.root_rank.below().unwrap();
        for root_child in &root_branch.children[0..rbc_len - 1] {
            dive_recursive(root_child, &self.payloads, false, String::new(), rank_below, f)?;
        }
        dive_recursive(&root_branch.children[rbc_len - 1], &self.payloads, true, String::new(), rank_below, f)
    }
}

/*
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

        eprintln!("{tree}");

        let expected = indoc! {"
            Root (Test)
            ├── test-1 0.075
            │   ├── test-11 0.15
            │   │   ├── test-111 0.1
            │   │   └── test-112 0.2
            │   └── test-12 0
            │       ├── test-121 0
            │       └── test-122 0
            └── test-2 0.9
                └── test-21 0.9
                    └── test-211 0.9
        "};
        assert_eq!(tree.to_string().as_str(), expected);
    }
}
*/
