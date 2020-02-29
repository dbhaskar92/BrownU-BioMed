#!/bin/bash

# Shell script to create an MPEG format movie from images 
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>

# Movie framerate
FPSEC="60"
MOVIENAME="my_movie"
SEPARATOR="_"

# Convert tif to png, if needed
# for f in *.tif; do convert "$f" "$(basename "$f" .tif).png"; done

CNT="0"
for f in *.png; do
	FRAME=$(echo $f | tail -c 9);
	NEWFILENAME=$MOVIENAME$SEPARATOR$FRAME
	mv "$f" "$NEWFILENAME";
	((CNT++)) 
done

LENGTH=$((CNT/FPSEC))
echo "Creating Movie - Frames: $CNT Estimated Length: $LENGTH seconds"

# Encode movie
INPUTIMAGES=$MOVIENAME$SEPARATOR
INPUTIMAGES+="%04d.png"
MOVIENAME+=".mp4"
avconv -framerate $FPSEC -i $INPUTIMAGES -b:v 8192k $MOVIENAME
