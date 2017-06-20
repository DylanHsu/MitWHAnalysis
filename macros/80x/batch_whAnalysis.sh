mkdir -p MitWHAnalsis/macros/80x
mv *.C MitWHAnalsis/macros/80x
mv *.h MitWHAnalsis/macros/80x
root -b -l -q MitWHAnalsis/macros/80x/whAnalysis.C+\(\"\",\"$3\",true,\"$1\",$2\)
