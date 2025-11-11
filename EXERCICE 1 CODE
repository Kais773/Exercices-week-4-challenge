% ====================================================================
% EXERCISE 1 - VARIANT B: THREE-PHASE STAR BALANCED LOAD
% ====================================================================

clear; clc; close all;

%% GIVEN DATA
UL = 230;           % Line voltage (V)
f = 60;             % Frequency (Hz)
R = 12;             % Resistance (Ω)
L = 0.05;           % Inductance (H)

%% CALCULATIONS
% Step 1: Inductive Reactance
XL = 2*pi*f*L;

% Step 2: Phase Impedance
Zph_mag = sqrt(R^2 + XL^2);
phi = atan2(XL, R);  % radians

% Step 3: Phase Voltage
Uph = UL/sqrt(3);

% Step 4: Phase Current
Iph = Uph/Zph_mag;

% Step 5: Phase angles
IA_angle = -rad2deg(phi);
IB_angle = -120 - rad2deg(phi);
IC_angle = 120 - rad2deg(phi);

% Step 6: Verify neutral current = 0
IA = Iph * exp(1j*deg2rad(IA_angle));
IB = Iph * exp(1j*deg2rad(IB_angle));
IC = Iph * exp(1j*deg2rad(IC_angle));
IN = IA + IB + IC;

% Step 7: Power calculations
cos_phi = cos(phi);
sin_phi = sin(phi);
Pph = Uph * Iph * cos_phi;
Qph = Uph * Iph * sin_phi;
Sph = Uph * Iph;
P3ph = 3 * Pph;
Q3ph = 3 * Qph;
S3ph = 3 * Sph;

% Step 8: Power factor correction
target_pf = 0.95;
if cos_phi < target_pf
    QC_total = P3ph * (tan(phi) - tan(acos(target_pf)));
    C_per_phase = QC_total / (3 * 2*pi*f * Uph^2);
else
    QC_total = 0;
    C_per_phase = 0;
end

%% PRINT RESULTS
fprintf('\n========== RESULTS ==========\n');
fprintf('XL = %.4f Ω\n', XL);
fprintf('|Zph| = %.4f Ω ∠ %.2f°\n', Zph_mag, rad2deg(phi));
fprintf('Uph = %.4f V\n', Uph);
fprintf('Iph = IL = %.4f A\n', Iph);
fprintf('|IN| = %.6f A (≈ 0)\n', abs(IN));
fprintf('\nPOWER PER PHASE:\n');
fprintf('Pph = %.2f W\n', Pph);
fprintf('Qph = %.2f VAR\n', Qph);
fprintf('Sph = %.2f VA\n', Sph);
fprintf('\nTOTAL 3-PHASE POWER:\n');
fprintf('P3φ = %.2f W = %.3f kW\n', P3ph, P3ph/1000);
fprintf('Q3φ = %.2f VAR = %.3f kVAR\n', Q3ph, Q3ph/1000);
fprintf('S3φ = %.2f VA = %.3f kVA\n', S3ph, S3ph/1000);
fprintf('Power Factor = %.4f = %.2f%%\n', cos_phi, cos_phi*100);
if cos_phi < target_pf
    fprintf('\nPF CORRECTION NEEDED:\n');
    fprintf('QC = %.2f kVAR\n', QC_total/1000);
    fprintf('C per phase = %.2f µF (star)\n', C_per_phase*1e6);
end
fprintf('=============================\n\n');

%% FIGURE 1: PHASOR DIAGRAM
figure('Position', [100, 100, 900, 700]);

% Voltage phasors
UAN = Uph * exp(1j*0);
UBN = Uph * exp(1j*deg2rad(-120));
UCN = Uph * exp(1j*deg2rad(120));

% Scale currents for visibility
scale = Uph/Iph * 0.8;
IA_scaled = IA * scale;
IB_scaled = IB * scale;
IC_scaled = IC * scale;

% Plot
hold on; grid on;
quiver(0, 0, real(UAN), imag(UAN), 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
quiver(0, 0, real(UBN), imag(UBN), 0, 'g', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
quiver(0, 0, real(UCN), imag(UCN), 0, 'b', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
quiver(0, 0, real(IA_scaled), imag(IA_scaled), 0, 'r--', 'LineWidth', 2);
quiver(0, 0, real(IB_scaled), imag(IB_scaled), 0, 'g--', 'LineWidth', 2);
quiver(0, 0, real(IC_scaled), imag(IC_scaled), 0, 'b--', 'LineWidth', 2);

text(real(UAN)*1.1, imag(UAN)*1.1, 'U_{AN}', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
text(real(UBN)*1.1, imag(UBN)*1.1, 'U_{BN}', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g');
text(real(UCN)*1.1, imag(UCN)*1.1, 'U_{CN}', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
text(real(IA_scaled)*1.1, imag(IA_scaled)*1.1, 'I_A', 'FontSize', 11, 'Color', 'r');
text(real(IB_scaled)*1.1, imag(IB_scaled)*1.1, 'I_B', 'FontSize', 11, 'Color', 'g');
text(real(IC_scaled)*1.1, imag(IC_scaled)*1.1, 'I_C', 'FontSize', 11, 'Color', 'b');

axis equal;
xlabel('Real', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imaginary', 'FontSize', 12, 'FontWeight', 'bold');
title('Phasor Diagram - Voltages & Currents', 'FontSize', 14, 'FontWeight', 'bold');
legend('U_{AN}', 'U_{BN}', 'U_{CN}', 'I_A', 'I_B', 'I_C', 'Location', 'best');
hold off;

%% FIGURE 2: POWER TRIANGLE
figure('Position', [150, 150, 900, 700]);
hold on; grid on;

plot([0 P3ph], [0 0], 'b', 'LineWidth', 3);
plot([P3ph P3ph], [0 Q3ph], 'r', 'LineWidth', 3);
plot([0 P3ph], [0 Q3ph], 'g', 'LineWidth', 3);
fill([0 P3ph P3ph 0], [0 0 Q3ph 0], [0.9 0.95 1], 'FaceAlpha', 0.3);

text(P3ph/2, -S3ph*0.08, sprintf('P = %.2f kW', P3ph/1000), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b', 'HorizontalAlignment', 'center');
text(P3ph*1.05, Q3ph/2, sprintf('Q = %.2f kVAR', Q3ph/1000), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
text(P3ph/2, Q3ph/2, sprintf('S = %.2f kVA', S3ph/1000), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g');

axis equal;
xlabel('Active Power P (W)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Reactive Power Q (VAR)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Power Triangle | PF = %.2f%%', cos_phi*100), 'FontSize', 14, 'FontWeight', 'bold');
legend('P (Active)', 'Q (Reactive)', 'S (Apparent)', 'Location', 'best');
hold off;

%% FIGURE 3: CIRCUIT DIAGRAM
figure('Position', [200, 200, 1000, 600]);
hold on; axis off; axis equal;

% Draw source
for i = 0:2
    y = 5 - i*2;
    rectangle('Position', [0.5 y-0.3 1 0.6], 'Curvature', 0.3, 'EdgeColor', 'k', 'LineWidth', 2);
    plot([0.7 1.3], [y y], 'k', 'LineWidth', 2);
    text(0.2, y, ['Phase ' char(65+i)], 'FontSize', 11, 'FontWeight', 'bold');
    plot([1.5 4.5], [y y], 'k', 'LineWidth', 2);
end

% Draw loads
for i = 0:2
    y = 5 - i*2;
    rectangle('Position', [5 y-0.25 1 0.5], 'EdgeColor', 'b', 'LineWidth', 2);
    text(5.5, y-0.5, 'Z_{ph}', 'FontSize', 10, 'HorizontalAlignment', 'center');
    plot([6 7], [y 3], 'k', 'LineWidth', 2);
end

plot(1, 3, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(7, 3, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot([1 7], [3 3], 'k--', 'LineWidth', 2);
text(4, 2.5, 'Neutral', 'FontSize', 11, 'FontWeight', 'bold');

text(5, 6.5, 'Three-Phase Star Connection', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(8.5, 3, sprintf('U_L = %.0f V\nU_{ph} = %.1f V\nI_{ph} = %.3f A', UL, Uph, Iph), 'FontSize', 11, 'BackgroundColor', [1 1 0.9]);

xlim([0 10]);
ylim([0 7]);
hold off;

%% FIGURE 4: PF CORRECTION (if needed)
if cos_phi < target_pf
    figure('Position', [250, 250, 1200, 500]);
    
    % Before correction
    subplot(1,2,1);
    hold on; grid on;
    plot([0 P3ph], [0 0], 'b', 'LineWidth', 3);
    plot([P3ph P3ph], [0 Q3ph], 'r', 'LineWidth', 3);
    plot([0 P3ph], [0 Q3ph], 'g', 'LineWidth', 3);
    fill([0 P3ph P3ph 0], [0 0 Q3ph 0], [1 0.9 0.9], 'FaceAlpha', 0.3);
    title(sprintf('BEFORE: PF = %.2f%%', cos_phi*100), 'FontSize', 13, 'Color', 'r');
    xlabel('P'); ylabel('Q');
    axis equal;
    hold off;
    
    % After correction
    subplot(1,2,2);
    hold on; grid on;
    Q_new = P3ph * tan(acos(target_pf));
    plot([0 P3ph], [0 0], 'b', 'LineWidth', 3);
    plot([P3ph P3ph], [0 Q_new], 'r', 'LineWidth', 3);
    plot([0 P3ph], [0 Q_new], 'g', 'LineWidth', 3);
    plot([P3ph P3ph], [Q3ph Q_new], 'm--', 'LineWidth', 3);
    fill([0 P3ph P3ph 0], [0 0 Q_new 0], [0.9 1 0.9], 'FaceAlpha', 0.3);
    text(P3ph*1.05, (Q3ph+Q_new)/2, sprintf('Q_C = %.2f kVAR', QC_total/1000), 'FontSize', 11, 'Color', 'm', 'FontWeight', 'bold');
    title(sprintf('AFTER: PF = %.2f%%', target_pf*100), 'FontSize', 13, 'Color', 'g');
    xlabel('P'); ylabel('Q');
    axis equal;
    hold off;
end

fprintf('All diagrams created successfully!\n');
