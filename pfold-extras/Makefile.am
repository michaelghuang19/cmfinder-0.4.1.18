DEFAULT_INCLUDES=-I.
LDADD=-lm
bin_PROGRAMS=addparen col2fasta col2newick col2psalign col2txtalign drawdot drawphyl extendstem findphyl mltree newick2col nohead scfg stdpair
addparen_SOURCES = addparen.c col.c file.c llist.c graph.c
col2fasta_SOURCES = col2fasta.c col.c file.c llist.c
col2newick_SOURCES = col2newick.c file.c llist.c phyl.c grammar.c optimize.c col.c align.c matrix.c edouble.c
col2psalign_SOURCES = col2psalign.c col.c file.c llist.c
col2txtalign_SOURCES = col2txtalign.c col.c file.c llist.c
drawphyl_SOURCES = drawphyl.c file.c llist.c phyl.c grammar.c optimize.c col.c align.c matrix.c edouble.c
extendstem_SOURCES = extendstem.c col.c file.c llist.c rna.c
findphyl_SOURCES = findphyl.c file.c llist.c phyl.c grammar.c optimize.c col.c align.c matrix.c edouble.c
newick2col_SOURCES = newick2col.c file.c llist.c phyl.c grammar.c optimize.c col.c align.c matrix.c edouble.c
nohead_SOURCES = nohead.c col.c file.c llist.c
stdpair_SOURCES = stdpair.c col.c file.c llist.c rna.c
drawdot_SOURCES = drawdot.c file.c llist.c
scfg_SOURCES = scfg.c phyl.c matrix.c col.c file.c llist.c grammar.c newcolprob.c edouble.c inout.c optimize.c align.c search.c
mltree_SOURCES = mltree.c file.c llist.c phyl.c grammar.c optimize.c col.c align.c matrix.c edouble.c

# don't build scfg_slow, since it requires compiling files with alternate flags -DETYPE=EDBL, e.g. to make matrix.o it's 'gcc -o matrix.o -c matrix.c', but to make matrix_slow.o it's 'gcc -DETYPE=EDBL -o matrix_slow.o -c matric.c'.  I don't need the scfg_slow command, and so I'd rather not figure out how to do this with automake
#scfg_slow_SOURCES = scfg.c phyl.c matrix_slow.c col.c file.c llist.c grammar_slow.c newcolprob_slow.c edouble_slow.c inout_slow.c optimize.c align.c search.c
