

function para = set_paras(mouse)
  
  para = struct;    %% maybe outside?
  para.f = 15;                 %% Hz frame rate
  para.nframe = 8989;
  
  
  if any(mouse==[884])
    para.totallength = 1600;     %% length of the linear track
    para.pos_offset = 1600;   %% offset of position values
    
    para.nSes = 15;
    
    %%% define timepoints of measurements
    para.t_s(1) = 0;     %% 09.25. AM
    para.t_s(2) = 4;     %% 09.25. PM
    para.t_s(3) = 24;     %% 09.26. AM
    para.t_s(4) = 28;     %% 09.26. PM
    para.t_s(5) = 48;     %% 09.27. AM
    para.t_s(6) = 52;     %% 09.27. PM
    para.t_s(7) = 120;     %% 09.30. AM
    para.t_s(8) = 124;     %% 09.30. PM
    para.t_s(9) = 148;     %% 10.01. PM
    para.t_s(10) = 168;     %% 10.02. AM
    para.t_s(11) = 172;     %% 10.02. PM
    para.t_s(12) = 192;     %% 10.03. AM
    para.t_s(13) = 196;     %% 10.03. PM
    para.t_s(14) = 216;     %% 10.04. AM
    para.t_s(15) = 220;     %% 10.04. PM
  
  elseif any(mouse==[245])
    %%% mouse 245   %% should the end and start position be contained in the "active" data? -> bias towards information during start/end position!?
    para.totallength = 600;     %% length of the linear track
    para.pos_offset = 300;   %% offset of position values
    
    para.nSes = 30;
    
    para.t_s(1) = 0;     %% 01.16. AM
    para.t_s(2) = 4;     %% 01.16. PM
    para.t_s(3) = 72;     %% 01.19. AM
    para.t_s(4) = 96;     %% 01.20. AM
    para.t_s(5) = 100;     %% 01.20. PM
    para.t_s(6) = 120;     %% 01.21. AM
    para.t_s(7) = 124;     %% 01.21. PM
    para.t_s(8) = 144;     %% 01.22. AM
    para.t_s(9) = 148;     %% 01.22. PM
    para.t_s(10) = 168;     %% 01.23. AM
    para.t_s(11) = 172;     %% 01.23. PM
    para.t_s(12) = 240;     %% 01.26. AM
    para.t_s(13) = 244;     %% 01.26. PM
    para.t_s(14) = 264;     %% 01.27. AM
    para.t_s(15) = 268;     %% 01.27. PM
    para.t_s(16) = 288;     %% 01.28. AM
    para.t_s(17) = 292;     %% 01.28. PM
    para.t_s(18) = 312;     %% 01.29. AM
    para.t_s(19) = 316;     %% 01.29. PM
    para.t_s(20) = 336;     %% 01.30. AM
    para.t_s(21) = 340;     %% 01.30. PM
    para.t_s(22) = 408;     %% 02.02. AM
    para.t_s(23) = 412;     %% 02.02. PM
    para.t_s(24) = 432;     %% 02.03. AM
    para.t_s(25) = 436;     %% 02.03. PM
    para.t_s(26) = 456;     %% 02.04. AM
    para.t_s(27) = 460;     %% 02.04. PM
    para.t_s(28) = 480;     %% 02.05. AM
    para.t_s(29) = 484;     %% 02.05. PM
    para.t_s(30) = 504;     %% 02.06. AM
  end
  
  
  para.nbin=80;                %% number of divisions of the linear track
  para.binwidth = para.totallength/para.nbin;
  
end