% include all sessions better that performed in the off case better thanchance
%% Session information
% hyper parameters were optmised on RSC data (VIP-Chr and VIP-Arch),
% thererfor there is a separate field called rmse_val which was hold out data that was used for the
% optimisation
sessionInfo = struct('type',{'VIP-Chr', 'VIP-Arch','VIP-Chr-HPC','VIP-Arch-HPC'},...% 'VIP-Chr-HPC','VIP-Arch-HPC'
                    'sessions',{{'m8058-20200530-00', 'm8059-20200530-00', 'm8059-20200623-00',...
                                'm8068-20211126-00', 'm8070-20211107-00', 'm8074-20211109-00',...
                                'm8043-20191213-00', 'm8058-20200527-01'},...% VIP-Chr
                                {'m8060-20200618-00', 'm8060-20200622-00', 'm8061-20200517-00',...
                                'm8063-20200618-00', 'm8063-20200624-00', 'm8071-20211107-00',...
                                'm8072-20211127-00', 'm8072-20211202-00', 'm8063-20200625-00',...
                                'm8061-20200523-00','m8061-20200519-00'},... % VIP-Arch
                                {'m8122-20230615-02','m8122-20230705-01','m8122-20230705-02',...
                                'm8123-20230615-01','m8123-20230616-01','m8123-20230704-01','m8123-20230705-02','m8123-20230705-01',...
                                'm8137-20240117-01','m8138-20240117-01','m8139-20240117-02','m8140-20240117-01'},...
                                {'m8117-20230705-02','m8118-20230705-02','m8120-20230705-02','m8121-20230705-02',...
                                 'm8083-20220407-00','m8098-20220709-00','m8117-20230615-01',...
                                 'm8118-20230615-01','m8118-20230705-01','m8119-20230615-01','m8119-20230616-01',...
                                 'm8120-20230615-01','m8120-20230705-01','m8121-20230615-01','m8121-20230616-01'}},...
                                  'mouse',{[4 5 5 1 3 2 6 4],[1 1 3 4 4 2 5 5 4 3 3],[1 1 1 2 2 2 2 2 3 4 5 6 ],[1 2 3 4 5 6 1 2 2 7 7 3 3 4 4]},...
                     'significance', {'right', 'left','right', 'left'});
                 
