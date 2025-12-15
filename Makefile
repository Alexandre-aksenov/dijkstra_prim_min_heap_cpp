target: 
	mkdir target

build: target
	g++ -o target/main.out main.cpp

run: target/main.out
	./target/main.out

clean:
	rm -f target/main.out
	