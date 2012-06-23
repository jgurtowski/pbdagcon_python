  $ cd $TESTDIR
  $ hcon.py r haplotype_test.fa htest_h1_result.fa --output ctest.fa
  $ blasr ctest_h1.fa htest_h1_result.fa 
  consensus_h1/0_3063 consensus_group_0_h1 0 0 -15315 100 0 3063 3063 0 3063 3063 64291
  $ blasr ctest_h2.fa htest_h2_result.fa 
  consensus_h2/0_3063 consensus_group_0_h2 0 0 -15315 100 0 3063 3063 0 3063 3063 64291
  $ blasr ctest_h1.fa ctest_h2.fa 
  consensus_h1/0_3063 consensus_h2 0 0 -15238 99.7715 0 3063 3063 0 3063 3063 64308
  $ rm ctest*
