inputfitslist := fitsList.in
path_fitsiolib := /opt/local/lib/
path_fitsioh := /opt/local/include/

fits2csv: fits2csv.cpp
	g++ -Wall -O2 -L$(path_fitsiolib) -I$(path_fitsioh) -lcfitsio fits2csv.cpp -o fits2csv

fits2csv_normalize: fits2csv_normalize.cpp
	g++ -Wall -O2 -L$(path_fitsiolib) -I$(path_fitsioh) -lcfitsio fits2csv_normalize.cpp -o fits2csv_normalize

exec:
	./fits2csv < $(inputfilename)

execnorm:
	./fits2csv_normalize < $(inputfilename)

clean:
	rm -f *.o readfits_csv

