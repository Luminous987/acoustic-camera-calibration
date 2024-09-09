function Mic_pos_err = compute_RMS_error(g)
Mic_pos_err = 0;

    for i = 2:g.M
        mic_idx = (i - 1) * 3;

        Mic_pos_err = Mic_pos_err + ((g.x_gt(mic_idx + 1) - g.x(mic_idx + 1))^2 + ...
        (g.x_gt(mic_idx + 2) - g.x(mic_idx + 2))^2 + (g.x_gt(mic_idx + 3) - g.x(mic_idx + 3))^2);
          
    end
    Mic_pos_err = (Mic_pos_err / g.M)^0.5;
end