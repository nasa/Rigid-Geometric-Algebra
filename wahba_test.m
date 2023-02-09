function wahba_test(opts)
%% Test Wahba
arguments
    opts.points (1,1) double = 1
    opts.point_sigma (1,1) double = 0
    opts.planes (1,1) double = 1
    opts.plane_sigma (1,1) = 0
    opts.directions (1,1) double = 1
    opts.direction_sigma (1,1) double = 0
    opts.lines (1,1) = 1
    opts.line_sigma (1,1) double = 0
    opts.motors (1,1) double = 1
    opts.motor_sigma (1,1) double = 0
    opts.rotation_sigma (1,1) double = 0
    opts.translation_sigma (1,1) double = 0
    opts.axis_sigma (1,1) double = 0
    opts.fullpose (1,1) logical = true
end
Qtru = unitize(rgamotor); Qtru.anti = true;
if ~opts.fullpose
    Qtru = weight(Qtru);
end
[~,~,vtru,mtru] = extract(Qtru);
nobs = opts.points + opts.planes + opts.directions + opts.lines + opts.motors;
for i = nobs:-1:1
    pherr = opts.rotation_sigma*randn;
    derr = opts.translation_sigma*randn;
    vper = vtru + opts.axis_sigma*randn(1,3);
    mper = mtru + opts.axis_sigma*randn(1,3);
    mper = reject(mper,vper);
    Lper = unitize(rgaline(vper,mper)); % perturbed axis (dir & mom)
    dR = unitize(rgamotor(pherr,0,Lper)); % rotation error
    dS = rgamotor(0,derr,Lper); % translation error
    dR.anti = true; dS.anti = true;
    dQ = rgamotor(dR*dS); % motor error
    Q = dQ*Qtru*~dQ;
    if opts.points
        M{i} = unitize(rgapoint); M{i}.anti = true;
        N{i} = rgapoint(Q*M{i}*~Q); N{i}.anti = true;
        if opts.point_sigma
            N{i}.m = N{i}.m + opts.point_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.points = opts.points - 1;
    elseif opts.planes
        M{i} = unitize(rgaplane); % anti by default
        N{i} = rgaplane(Q*M{i}*~Q);
        if opts.plane_sigma
            N{i}.m = N{i}.m + opts.plane_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.planes = opts.planes - 1;
    elseif opts.directions
        M{i} = rgaline(unit(randn(1,3)),zeros(1,3)); M{i}.anti = true;
        N{i} = rgaline(Q*M{i}*~Q); N{i}.anti = true;
        if opts.direction_sigma
            N{i}.m = N{i}.m + opts.direction_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.directions = opts.directions - 1;
    elseif opts.lines
        M{i} = unitize(rgaline); M{i}.anti = true;
        N{i} = rgaline(Q*M{i}*~Q); N{i}.anti = true;
        if opts.line_sigma
            N{i}.m = N{i}.m + opts.line_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.lines = opts.lines - 1;
    elseif opts.motors
        M{i} = unitize(rgamotor); M{i}.anti = true;
        N{i} = rgamotor(Q*M{i}*~Q); N{i}.anti = true;
        if opts.motor_sigma
            N{i}.m = N{i}.m + opts.motor_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.motors = opts.motors - 1;
    end
end

Qest = wahba(M,N);
Qest.anti = true;
disp('Qtru = ')
disp(Qtru)
disp('Qest =')
disp(Qest)
disp('Qest - Q =')
disp(Qest - Q)
[phit,dt,vt,mt] = extract(Qtru);
[phie,de,ve,me] = extract(Qest);
disp('Rotation error [deg] (phie - phit) =')
disp(180/pi*(phie-phit))
disp('Translation error (de - dt) =')
disp(de - dt)
disp('True Direction & Moment of Screw Axis = ')
disp([vt mt])
disp('Est Direction & Moment of Screw Axis = ')
disp([ve me])

end

function a2 = reject(a,b)
a2 = a - dot(a,b)/dot(b,b)*b;
end