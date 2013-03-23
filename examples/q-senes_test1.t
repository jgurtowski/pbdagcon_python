  $ cd $TESTDIR
  $ q-sense.py d consensus_test.fa --output ctest.fa
  $ blasr ctest.fa consensus_result.fa 
  consensus/0_8198 consensus 0 0 -40625 99.9754 0 8129 8129 71 8198 8198 170641
  $ blasr ctest.fa consensus_reference.fa 
  consensus/0_8198 ref 0 0 -40980 99.9756 898 9098 9098 0 8198 8198 172132
  $ rm ctest*
