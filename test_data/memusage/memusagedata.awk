#!/bin/awk -f

BEGIN{
	FS = "\t"
	OFS = "\t"
}
{
	if ($0 ~ /^Max Mem used/)
	{
		gsub(/ /, "", $0)
		split($0, a, ":")
		split(a[2], b, "(")
		if (substr(b[1], length(b[1]), 1) == "M")
		{
			print substr(b[1], 1, length(b[1]) - 1) / 1000
		}
		else
		{
			print substr(b[1], 1, length(b[1]) - 1)
		}
	}
}