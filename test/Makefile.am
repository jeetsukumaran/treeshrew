AM_CPPFLAGS = -std=c++11 @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

# AUTOMAKE_OPTIONS = subdir-objects

treeshrewtestsdir = $(prefix)/opt/@PACKAGE@/test
# treeshrewtestsdir = $(builddir)/bin

dist_treeshrewtests_SCRIPTS = \
	scripts/run-tests.py

nobase_dist_treeshrewtests_DATA = \
	data/basic/bird_orders.nex \
	data/basic/d20121031.txt \
	data/basic/primates.beast.mcct.medianh.newick.tre \
	data/basic/primates.beast.mcct.medianh.tre \
	data/basic/primates.chars.fasta \
	data/basic/primates.chars.nexus \
	data/basic/pythonidae.beast.summary.tre \
	data/basic/pythonidae.chars.fasta \
	data/basic/pythonidae.chars.nexus \
	data/basic/pythonidae.postorder.newick \
	data/basic/pythonidae.preorder.newick \
	data/basic/pythonidae.tree.newick

COMMON_TREE_SRC = \
	../src/utility.hpp \
	../src/utility.cpp \
	../src/tree.hh \
	../src/genetree.hpp \
	../src/genetree.cpp \
	../src/dataio.hpp

COMMON_TEST_SRC = \
	src/testutils.hpp \
	src/testutils.cpp

treeshrewtests_PROGRAMS = \
	phylogenetic_likelihood_calc \
	read_tree \
	tree_postorder_iter \
	tree_preorder_iter \
	tree_leaf_iter \
	tree_internal_iter \
	tree_child_iter \
	read_dna_alignment \
	score_phylogenetic_tree

phylogenetic_likelihood_calc_SOURCES = \
	src/phylogenetic_likelihood_calc.cpp

read_tree_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/read_tree.cpp

tree_postorder_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_postorder_iter.cpp

tree_preorder_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_preorder_iter.cpp

tree_leaf_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_leaf_iter.cpp

tree_internal_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_internal_iter.cpp

tree_child_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_child_iter.cpp

read_dna_alignment_SOURCES = \
	$(COMMON_TEST_SRC) \
	../src/utility.hpp \
	../src/utility.cpp \
	../src/character.hpp \
	../src/character.cpp \
	../src/tree.hh \
	../src/genetree.hpp \
	../src/genetree.cpp \
	../src/dataio.hpp \
	src/read_dna_alignment.cpp

score_phylogenetic_tree_SOURCES = \
	$(COMMON_TEST_SRC) \
	../src/utility.hpp \
	../src/utility.cpp \
	../src/character.hpp \
	../src/character.cpp \
	../src/tree.hh \
	../src/genetree.hpp \
	../src/genetree.cpp \
	../src/dataio.hpp \
	src/score_phylogenetic_tree.cpp

