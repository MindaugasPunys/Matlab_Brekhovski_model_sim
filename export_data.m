function export_data(data, var_name)
% Check if data is complex
is_complex = ~isreal(data);

if is_complex
    fid_re = fopen([var_name, '_re.txt'], 'w');
    fid_im = fopen([var_name, '_im.txt'], 'w');

    fprintf(fid_re, '%f\n', real(data));
    fprintf(fid_im, '%f\n', imag(data));

    fclose(fid_re);
    fclose(fid_im);
else
    fid = fopen([var_name, '.txt'], 'w');
    fprintf(fid, '%f\n', data);
    fclose(fid);
end

end