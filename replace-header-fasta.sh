#! /bin/zsh

	for file in 'ls "$1"'
		do
		sed 's|^>.*REF|>"$file|' <$file >$file
		done
