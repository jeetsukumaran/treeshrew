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
	../src/character.hpp \
	../src/character.cpp \
	../src/tree.hpp \
	../src/genetree.hpp \
	../src/genetree.cpp \
	../src/statespace.hpp \
	../src/statespace.cpp \
	../src/dataio.hpp

COMMON_TEST_SRC = \
	src/testutils.hpp \
	src/testutils.cpp

treeshrewtests_PROGRAMS = \
	check_gsl_installation \
	phylogenetic_likelihood_calc \
	read_tree \
	tree_postorder_iter \
	tree_preorder_iter \
	tree_leaf_iter \
	tree_internal_iter \
	tree_child_iter \
	read_dna_sequences \
	score_phylogenetic_tree \
	benchmark_phylogenetic_tree \
	calc_hamming_distance

check_gsl_installation_SOURCES = \
	src/check_gsl_installation.cpp

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

read_dna_sequences_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/read_dna_sequences.cpp

score_phylogenetic_tree_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/score_phylogenetic_tree.cpp

benchmark_phylogenetic_tree_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/benchmark_phylogenetic_tree.cpp

calc_hamming_distance_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/calc_hamming_distance.cpp

