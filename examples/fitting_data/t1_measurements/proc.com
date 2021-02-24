#/usr/bin/csh

foreach dir (`ls -d data/*.fid`)
    echo "Processing files in $dir"
    cd $dir
    cp ../../xy_s3e.com .
    ./xy_s3e.com
    cd ../..
end
