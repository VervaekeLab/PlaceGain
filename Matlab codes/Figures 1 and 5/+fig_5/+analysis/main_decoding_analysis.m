% load_session_info()

load_dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/sData/VIP-Chr-HPC/';
save_dirct ='/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/15secOpto/Decoding/VIP-Chr-HPC';

for p = 1:length(sessionInfo(3).sessions)
    file = sessionInfo(3).sessions{p};

    fig_5.analysis.decode_dff()
end

