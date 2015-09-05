all: animation chimax large-L-one-gauss large-L-one-Mchi l-one

chimax: chimax.cc ~/include/characters.h ~/include/slint.h
	g++ -O3 -o chimax chimax.cc -pthread -lfftw3 -std=c++11
l-one: l-one.cc ~/include/characters.h ~/include/slint.h
	g++ -O3 -o l-one l-one.cc -std=c++0x -pthread -lfftw3 -std=c++11
large-L-one-Mchi: large-L-one-Mchi.cc ~/include/characters.h ~/include/slint.h
	g++ -O3 -o large-L-one-Mchi large-L-one-Mchi.cc -std=c++11 -pthread -lfftw3 -std=c++11
large-L-one-gauss: large-L-one-gauss.cc ~/include/characters.h ~/include/slint.h
	g++ -O3 -o large-L-one-gauss large-L-one-gauss.cc -std=c++11 -pthread -lfftw3 -std=c++11
animation: animation.cc ~/include/characters.h ~/include/slint.h
	g++ -std=c++11 -O3 -o animation animation.cc -pthread -lfftw3 -lSDL2 -std=c++11

clean:
	-rm chimax
	-rm l-one
	-rm large-L-one-Mchi
	-rm large-L-one-gauss
	-rm animation
