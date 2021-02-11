
# TODO: change the libPath to your specific path
libPath="/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/lib"

cmakePath="$libPath/cmake-3.17.2"
gccPath="$libPath/gcc-9.1.0"
CLI11Path="$libPath/CLI11-1.9.0"
spdlogPath="$libPath/spdlog-1.5.0"
libdeflatePath="$libPath/libdeflate-1.5"
htslibPath="$libPath/htslib-1.9"
yggPath="$libPath/ygg-master"
doctestPath="$libPath/doctest-2.3.7"
fftwPath="$libPath/fftw-3.3.8"

binPath="/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin"
export PATH="$gccPath/bin:$cmakePath/bin:$binPath"
#echo $PATH

export LD_LIBRARY_PATH="$fftwPath/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$libdeflatePath/lib:$htslibPath/lib:$gccPath/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/lib64"
export LIBRARY_PATH=$LD_LIBRARY_PATH
#echo $LD_LIBRARY_PATH

export C_INCLUDE_PATH="$fftwPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$doctestPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$libdeflatePath/include:$yggPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$htslibPath/include:$gccPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$CLI11Path/include:$spdlogPath/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH
#echo $C_INCLUDE_PATH

export CC="$gccPath/bin/gcc"
export CXX="$gccPath/bin/g++"

absPath(){
    relativePath=$1
    mkdir -p $relativePath
    if [ ${relativePath:0:1} == "/" ]
    then
        echo $relativePath
        return 0
    fi
    ap="$(cd $relativePath; pwd)"
    echo "$ap"
}

#srcPath=/hwfssz1/ST_BIGDATA/USER/zhaofuxiang/git/handle_bam
srcPath="$(cd $(dirname $(dirname $0)); pwd)"
buildPath="$(mktemp -d -p /dev/shm $(basename $srcPath).XXXXXXXXXX)"
#echo $buildPath
buildPath="$(absPath $buildPath)"
if [ -e "$buildPath" ]
then
    rm -rf $buildPath
fi
mkdir -p $buildPath

cd $buildPath
#echo $buildPath

installPath="$buildPath/install"
mkdir -p $installPath/bin $installPath/lib

timeStart=$(date +%s)
#cmake -DCMAKE_CXX_COMPILER=$gccPath/bin/c++ $srcPath
#make clean
test="OFF"
if [ $# == 1 ]
then
    if [ $1 == "ON" ]
    then
        test="ON"
    fi
fi

cmake $srcPath -DINSTALL_PATH=$installPath -DUNITTEST=$test

thread=$(grep -c ^processor /proc/cpuinfo)
make -j $thread #VERBOSE=1


#cp $buildPath/bin/handleBam $installPath/bin
#cp $buildPath/bin/handleBam /hwfssz5/ST_BIGDATA/USER/zhaofuxiang/test_handleBam/handleBam


install(){
    libFile="$1"
    tlib="$installPath/lib"
    if [ -d "$tlib" ]
    then
        if [ -e "$libFile" ]
        then
            echo "Installing: $libFile into $tlib"
            cp $libFile $tlib
        else
            echo "[WARN] File not exists: $libFile"
        fi
    fi
}

install $htslibPath/lib/libhts.so.2
install $libdeflatePath/lib/libdeflate.so.0

cd $srcPath
cp -R $installPath $srcPath
if [ -e "$buildPath" ]
then
    rm -rf $buildPath
fi

timeEnd=$(date +%s)
secs=$(($timeEnd - $timeStart))
mins=$(($secs/60))
secs=$(($secs%60))
echo "Cost: $mins Mins $secs Secs"
