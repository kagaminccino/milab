#! /bin/bash

python add_time.py

SRC=`find . -name "*.png" | sort -n`
i=0
for s in $SRC; do
    dst=$(printf ./cap%05d.png $i)
    if [ "$s" != "$dst" ]; then
        mv $s $(printf cap%05d.png $i)
    fi
    i=$(($i + 1))
done
# ffmpeg -framerate 50 -i cap%05d.png -c:v libx264 -pix_fmt yuv420p video.mp4
convert -delay 2 -loop 0 cap*.png video.gif
rm -f *.xwd
rm -f *.png
