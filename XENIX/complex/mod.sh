:
for i in *.[ch]; do
	cp "${i}" "${i}.O"
	eval `awk 'BEGIN { printf("sed "); } { printf("-e '\''s/%s/%s/g'\'' ", $1, $2); }' <symbols` <"${i}" >"${i}.N"
	mv "${i}.N" "${i}"
done
