  $ cd $TESTDIR
  $ gcon.py r consensus_test.fa consensus_result.fa --output ctest.fa
  $ blasr ctest.fa consensus_result.fa 
  consensus/0_8129 consensus 0 0 -40645 100 0 8129 8129 0 8129 8129 170677
  $ blasr ctest.fa consensus_reference.fa 
  consensus/0_8129 ref 0 0 -40645 100 969 9098 9098 0 8129 8129 170677
  $ rm ctest*
