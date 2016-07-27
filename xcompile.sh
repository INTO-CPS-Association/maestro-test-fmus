#!/bin/bash
set -e

threads=4

function compileDarwin64
{

		echo Building Darwin 64 .dylib
		B=$1/build/darwin64

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/osx-gcc.cmake` -DOSXCROSS_ROOT=$OSXCROSS_ROOT

		make -C $B -j$threads

}

function compileWin32
{

		echo Building Win32 .dll
		B=$1/build/win32

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/cmake-toolchains/Toolchain-Ubuntu-mingw32.cmake`

		make -C $B -j$threads

}

function compileWin64
{

		echo Building Win64 .dll
		B=$1/build/win64

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/cmake-toolchains/Toolchain-Ubuntu-mingw64.cmake`

		make -C $B -j$threads

}

function compileLinux64
{
		echo Building Linux x64 .so
		B=$1/build/linux64
		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1

		make -C $B -j$threads
}

function compileLinux32
{
		echo Building Linux x32 .so
		B=$1/build/linux32

		rm -rf $B
		mkdir -p $B

		cmake -B$B -H$1 -DLINK_FLAGS="-m32" -DCFLAGS="-m32"

		make -C $B -j$threads
}


function overtureToolWrapper
{
		wget http://overture.au.dk/into-cps/vdm-tool-wrapper/development/latest/vdm-tool-wrapper.zip -O vdm-tool-wrapper.zip

		unzip -l vdm-tool-wrapper.zip


    mkdir -p $1/sources/
    touch $1/sources/none.txt
		cp $1/modelDescription.xml $1/resources/modelDescription.xml
		
		# resources
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/resources/config.txt -d $1/resources/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/resources/fmi-interpreter-*-jar-with-dependencies.jar -d $1/resources/

		bindir=$1/build
		mkdir -p $bindir/{darwin64,win32,win64,linux32,linux64}
		#tool wrapper binaries
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/darwin64/vdm-tool-wrapper.dylib -d $bindir/darwin64/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/win32/vdm-tool-wrapper.dll -d $bindir/win32/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/win64/vdm-tool-wrapper.dll -d $bindir/win64/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/linux32/vdm-tool-wrapper.so -d $bindir/linux32/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/linux64/vdm-tool-wrapper.so -d $bindir/linux64/
		#rename
		
		
		mv $bindir/darwin64/*.dylib $bindir/darwin64/$2.dylib
		mv $bindir/linux64/*.so $bindir/linux64/$2.so
		mv $bindir/linux32/*.so $bindir/linux32/$2.so
		mv $bindir/win64/*.dll $bindir/win64/$2.dll
		mv $bindir/win32/*.dll $bindir/win32/$2.dll
}

function assemble
{

		B=$1/build/

		mkdir -p $B/fmu/{binaries,resources,sources}

		mkdir -p $B/fmu/binaries/{darwin64,win32,win64,linux32,linux64}

		echo Copying files...
		cp $1/modelDescription.xml $B/fmu/

		cp $1/sources/*.* $B/fmu/sources/
		cp $1/resources/*.* $B/fmu/resources/ 2>/dev/null || :

		BIN=$B/fmu/binaries
		cp $B/darwin64/*.dylib $BIN/darwin64/
		cp $B/linux64/*.so $BIN/linux64/
		cp $B/linux32/*.so $BIN/linux32/
		cp $B/win64/*.dll $BIN/win64/
		cp $B/win32/*.dll $BIN/win32/

		echo Zipping...

		curdir="$PWD"

		cd $B/fmu/
		zip -r ../$name.fmu .

		cd $curdir
}

function xcompile
{
		D=$1
		echo Compiling using CMake and make
		compileDarwin64 $D
		compileLinux64 $D
		compileLinux32 $D
		compileWin64 $D
		compileWin32 $D
}


for D in `find . -maxdepth 1 -type d ! -path . ! -path *.git ! -path */includes ! -path */toolchains`
do
		echo "Folder is ${D}"
		name=`echo $D| sed 's|\./||g'`
		echo "Cleadning old builds ..."
		rm -rf $D/build/
		if [ -e "$D/mode.txt" ] 
		then
				if grep -q OVERTURE_TOOL_WRAPPER "$D/mode.txt"; then
						echo "Overture tool wrapper"
						overtureToolWrapper $D $name
				else
						xcompile $D
				fi
		else
				xcompile $D
		fi
		assemble $D $name
done

