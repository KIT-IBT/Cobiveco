function exportResult(o)

if ~o.available.apicobasal
    o.computeApicobasal;
end
o.printStatus('Exporting result...');
t = toc;

if o.cfg.exportLevel > 0
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));

    R = o.R;
    save(sprintf('%sR.mat', o.cfg.outPrefix), 'R');

    if o.cfg.exportLevel > 2
        result_R = o.result;
        result_R.points = [result_R.points ones(size(result_R.points,1),1)]*R';
        result_R.points(:,end) = [];
        vtkWrite(result_R, sprintf('%sresult_R.vtu', o.cfg.outPrefix));
    end
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);

end