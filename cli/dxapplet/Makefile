all:
	cd ../.. && cmake -DCMAKE_BUILD_TYPE=Profile -Dtest=ON && make -j4 && ./unit_tests
	cp ../../glnexus_cli resources/usr/local/bin/

.PHONY: all
