function Analysis = ParameterSweepAnalysis(X, Y, K, Lambda, P)

    
    %%
    [Lambda, P] = meshgrid(Lambda, P);
    Lambda = Lambda(:); 
    P = P(:);
    
    %%
    
    Temp_dat = zeros(length(Lambda), length(X));
    Velo_dat = zeros(length(Lambda), length(X));
    Accel_dat = zeros(length(Lambda), length(X));
    
    for k = 1:length(Lambda)
    
        [~, ~, Temp_f, ~, ~] = NonStandardFourierSeries(Y, X, P(k), K, Lambda(k));
        Temp = Temp_f(X);
        Velo = interp1(X(2:end-1), (Temp(3:end) - Temp(1:end-2))./(X(3:end) - X(1:end-2)), X);
        Accel = interp1(X(2:end-1), (Velo(3:end) - Velo(1:end-2))./(X(3:end) - X(1:end-2)), X);
    
        Temp_dat(k, :) = Temp;
        Velo_dat(k, :) = Velo;
        Accel_dat(k, :) = Accel;
        disp(round(k/length(Lambda)*100, 2))
    
    end
    
    
    
    Analysis.Temp_dat = Temp_dat;
    Analysis.Velo_dat = Velo_dat;
    Analysis.Accel_dat = Accel_dat;
    
    Analysis.Temp_avg = mean(Temp_dat);
    Analysis.Velo_avg = mean(Velo_dat);
    Analysis.Accel_avg = mean(Accel_dat);
    
    Analysis.Temp_std = std(Temp_dat);
    Analysis.Velo_std = std(Velo_dat);
    Analysis.Accel_std = std(Accel_dat);
    
    disp('Data Collected')
    
    %%
    
    
    [X_space, Xi_space, pdf_matrix] = ComputeDensity(X(2:end-1), ...
                                                     Velo_dat(:, 2:end-1), ...
                                                     [-0.01, 0.04], 600);
    
    Analysis.VeloPDF1D.X_space = X_space;
    Analysis.VeloPDF1D.Xi_space = Xi_space; 
    Analysis.VeloPDF1D.pdf_matrix = pdf_matrix;
    
    disp('VeloPDF1 finishd')
    
    %%
    
    [X_space, Xi_space, pdf_matrix] = ComputeDensity(X(3:end-2), ...
                                                     Accel_dat(:, 3:end-2), ...
                                                     [-0.0015, 0.0035], 600);
    
    Analysis.AccelPDF1D.X_space = X_space;
    Analysis.AccelPDF1D.Xi_space = Xi_space; 
    Analysis.AccelPDF1D.pdf_matrix = pdf_matrix;
    
    disp('AccelPDF1 finishd')
    
    %%
    [Tempgrid, Velogrid, pdfTV] = ComputePhaseSpaceDensity(Temp_dat(:, 2:end-1), Velo_dat(:, 2:end-1), 400);
    Analysis.TVpdf.TempGrid = Tempgrid; 
    Analysis.TVpdf.VeloGrid = Velogrid; 
    Analysis.TVpdf.PDF = pdfTV;
    disp('2D Temp-Velo PDF finished')
    
    %%
    [Velogrid, Accelgrid, pdfVA] = ComputePhaseSpaceDensity(Velo_dat(:, 3:end-2), Accel_dat(:, 3:end-2), 400);
    Analysis.VApdf.VeloGrid = Velogrid; 
    Analysis.VApdf.AccelGrid = Accelgrid; 
    Analysis.VApdf.PDF = pdfVA;
    disp('2D Velo-Accel PDF finished')



end