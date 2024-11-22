rem %1 - $(SolutionDir)
rem %2 - $(Configuration)

cd %1\3rd_party

@echo "nasm"

if exist nasm/nasm.exe (
	@echo "nasm.exe already exists"
	cd nasm
) else (
	rmdir /S /Q nasm
	mkdir nasm
	cd nasm
	curl -L --ssl-no-revoke https://github.com/refresh-bio-dependencies/nasm/releases/download/v2.16.01/nasm-2.16.01-win64.zip --output nasm-2.16.01-win64.zip
	tar -xf nasm-2.16.01-win64.zip --strip-components 1
)

set PATH=%PATH%;%cd%
cd ..

@echo %PATH%

@echo "zlib-ng"
cd zlib-ng 
cmake -B build-vs -S . -DZLIB_COMPAT=ON 
cmake --build build-vs --config %2
cd ..


@echo "isa-l"
cd isa-l
nmake -f Makefile.nmake
