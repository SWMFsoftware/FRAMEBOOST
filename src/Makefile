test_boosted_frame.exe: 
	gfortran -o test_boosted_frame.exe \
		ModUtility.f90 ModExactRS.f90 ModBoostedFrame.f90 ModPlotFile.f90  \
		test_boosted_frame.f90

test: test_boosted_frame.exe
	./test_boosted_frame.exe

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f *.out *.exe
