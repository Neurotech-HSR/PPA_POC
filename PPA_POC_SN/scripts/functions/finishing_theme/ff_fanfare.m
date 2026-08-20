function ff_fanfare()
% Hope is not copyrighted
    fs  = 44100;
    bpm = 150;
    q   = 60/bpm;
    tri = q/3;   % terzina (rinominato da t per evitare conflitto)
    e   = q/2;   % ottavo
    s   = q/4;   % sedicesimo (non usato ora, lasciato per reference)
    h   = q*2;
    f   = q/3;

    E5 = 659.25;  B5 = 987.77;
    C5 = 523.25;  G5 = 783.99;
    D5 = 587.33;  A5 = 880.00;

    % [root, fifth, durata]  (0 = pausa)
    seq = [
        E5, B5, tri;
        E5, B5, tri;
        E5, B5, tri;
        E5, B5, q;
        C5, G5, q;
        D5, A5, q;
        E5, B5, f;
        0,  0,  e;    % pausa ottavo
        D5, A5, e;
        E5, B5, h;
    ];

    audio = [];
    for i = 1:size(seq,1)
        root = seq(i,1);
        fift = seq(i,2);
        dur  = seq(i,3);
        n    = round(dur * fs);
        tvec = (0:n-1) / fs;

        if root == 0
            seg = zeros(1, n);
        else
            seg = 0.35*sin(2*pi*root*tvec) + 0.25*sin(2*pi*fift*tvec);
            att = min(round(0.008*fs), n);
            rel = min(round(0.05*fs),  n);
            env = ones(1,n);
            env(1:att) = linspace(0,1,att);
            env(end-rel+1:end) = linspace(1,0,rel);
            seg = seg .* env;
        end

        audio = [audio, seg, zeros(1,round(0.01*fs))]; %#ok<AGROW>
    end
    sound(audio, fs);
end