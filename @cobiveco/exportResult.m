function exportResult(o)

% Exports the results.
%
% exportResult(o)
%
% Input:
%   o: instance of the class cobiveco
%
% Output:
%   o.R, object of class cobiveco [double, size 4x4] (For details see cobiveco class documentation)
%   o.leftRightAx, object of class cobiveco [double, size 1x3] (For details see cobiveco class documentation)
%   o.antPostAx, object of class cobiveco [double, size 1x3] (For details see cobiveco class documentation)
%   o.longAx, object of class cobiveco [double, size 1x3] (For details see cobiveco class documentation)

o.printStatus('Exporting result...');
t = toc;

if o.cfg.exportLevel > 0
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));

    R = o.R;
    save(sprintf('%sR.mat', o.cfg.outPrefix), 'R');

    if o.cfg.exportLevel > 2
        leftRightAx = o.leftRightAx;
        save(sprintf('%sleftRightAx.mat', o.cfg.outPrefix), 'leftRightAx');
        antPostAx = o.antPostAx;
        save(sprintf('%santPostAx.mat', o.cfg.outPrefix), 'antPostAx');
        longAx = o.longAx;
        save(sprintf('%slongAx.mat', o.cfg.outPrefix), 'longAx');
        resultR = o.result;
        resultR.points = [resultR.points ones(size(resultR.points,1),1)]*R';
        resultR.points(:,end) = [];
        vtkWrite(resultR, sprintf('%sresultR.vtu', o.cfg.outPrefix));
    end
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);

end

