#!/bin/bash
set -e

function compileDarwin
{

		echo Building Darwin .dylib
		B=$1/build/darwin

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/osx-gcc.cmake` -DOSXCROSS_ROOT=$OSXCROSS_ROOT

		make -C $B

}

function compileWin32
{

		echo Building Win32 .dll
		B=$1/build/win32

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/cmake-toolchains/Toolchain-Ubuntu-mingw32.cmake`

		make -C $B

}

function compileWin64
{

		echo Building Win64 .dll
		B=$1/build/win64

		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1 -DCMAKE_TOOLCHAIN_FILE=`readlink -f toolchains/cmake-toolchains/Toolchain-Ubuntu-mingw64.cmake`

		make -C $B

}

function compileLinux64
{
		echo Building Linux x64 .so
		B=$1/build/linux64
		rm -rf $B
		mkdir -p $B

		cmake  -B$B -H$1

		make -C $B
}

function compileLinux32
{
		echo Building Linux x32 .so
		B=$1/build/linux32

		rm -rf $B
		mkdir -p $B

		cmake -B$B -H$1 -DLINK_FLAGS="-m32" -DCFLAGS="-m32"

		make -C $B		
}


function overtureToolWrapper
{
		wget http://overture.au.dk/into-cps/vdm-tool-wrapper/development/latest/vdm-tool-wrapper.zip -O vdm-tool-wrapper.zip


		# resources
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/resources/config.txt -d $1/resources/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/resources/fmi-interpreter-0.0.1-SNAPSHOT-jar-with-dependencies.jar -d $1/resources/

		mkdir -p $1/binaries/
		#tool wrapper binaries
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/darwin64/libshmfmu.dylib -d $1/binaries/darwin64/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/win32/libshmfmu.dll -d $1/binaries/win32/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/win64/libshmfmu.dll -d $1/binaries/win64/
		unzip -o -j vdm-tool-wrapper.zip vdm-tool-wrapper/binaries/linux64/libshmfmu.so -d $1/binaries/linux64/

}

function assemble
{

		B=$1/build/

		mkdir -p $B/fmu/{binaries,resources,sources}

		mkdir -p $B/fmu/binaries/{darwin64,win32,win64,linux32,linux64}

		echo Copying files...
		cp $1/modelDescription.xml $B/fmu/

		cp $1/sources/*.* $B/fmu/sources/

		BIN=$B/fmu/binaries
		cp $B/darwin/*.dylib $BIN/darwin64/
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
		compileDarwin $D
		compileLinux64 $D
		compileLinux32 $D
		compileWin64 $D
		compileWin32 $D
}


for D in `find . -maxdepth 1 -type d ! -path . ! -path *.git ! -path */includes ! -path */toolchains`
do
		echo "Folder is ${D}"
		name=`echo $D| sed 's|\./||g'`
		echo $D
		if [ -e "$D/mode.txt" ] 
		then
				if grep -q OVERTURE_TOOL_WRAPPER "$D/mode.txt"; then
						echo "Overture tool wrapper"
						overtureToolWrapper $D
				else
						xcompile $D
				fi
		else
				xcompile $D
		fi
		assemble $D $name
done

