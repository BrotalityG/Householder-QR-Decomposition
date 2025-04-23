% Branden Stahl
% Ma 432 Final Project

clc;
clear;
close all;

% Matrix A
A = [1 1 0; 1 0 1; 0 1 1]; % This is our example, same as worked in class

% Our method using householder transformations
[Q, R] = hqr(A);

% what it should be (signs are insignificant for householder)
[Q1, R1] = qr(A);

fprintf("Note that for QR decompostion, signs are insignificant except with respect to Q & R.\n")

%% Is it correct?
QR = Q*R;

% Mask numbers with limits of themselves, due to computational rounding errors
epsilon = 10^-6;
QR = round(QR*(1/epsilon))*epsilon;  % Apply

correct = isequal(A, QR);
if correct == 1
    fprintf("The QR decomposition is correct.\nQ:\n");
    disp(round(Q*(1/epsilon))*epsilon);
    fprintf("R:\n")
    disp(round(R*(1/epsilon))*epsilon);
    fprintf("\nA:\n")
    disp(A);
    fprintf("QR:\n")
    disp(QR);
else
    fprintf("The QR decomposition is incorrect.\n");
end

%% Functions
function [Q, R] = hqr(A) % Our function to do QR decomposition
    [n,p] = size(A);
    m = min(n,p);

    Q = eye(m);
    R = A;

    for k = 1:n % QR decomp with householder is iterative
        x = R(k:n, k);
        e1 = zeros(length(x), 1);
        e1(1) = norm(x);
    
        % solving for householder matrix
        v = x - sign(x(1)) * e1;

        nv = norm(v);
        if nv ~= 0
            v = v / nv;
        else
            v = 0;
        end

        Hk = eye(n-k+1) - 2 * (v * v');
        H = eye(n);
        H(k:n, k:n) = Hk;

        % applying householder transformation
        R = H*R;
        Q = Q*H';
    end

    % Testing properties of householder: Hermitian, Unitary,
    % Involutory, and determinant is equal to +1 or -1 due to sign
    % issues within MATLAB.
    if isequal(H, H') & isequal(inv(H), H') & isequal(H, inv(H)) & abs(det(H)) == 1
        fprintf('A valid householder matrix was produced.\n');
    end
end
