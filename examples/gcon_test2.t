  $ cd $TESTDIR
  $ gcon.py r consensus_test.fa consensus_reference.fa --output ctest.fa
  $ blasr ctest.fa consensus_result.fa 
  consensus/0_8532 consensus 0 0 -40640 100 0 8128 8129 404 8532 8532 170656
  $ blasr ctest.fa consensus_reference.fa 
  consensus/0_8532 ref 0 0 -42650 99.9883 566 9097 9098 0 8532 8532 179144
  $ rm ctest*
