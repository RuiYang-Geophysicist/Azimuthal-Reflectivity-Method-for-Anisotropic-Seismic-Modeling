%% compile_mex.m
% Compile AzRM C++ code into MATLAB MEX function
%
% Usage:
%   compile_mex          % Auto-detect settings
%   compile_mex('debug') % Debug build
%
% Prerequisites:
%   - MATLAB with mex command configured
%   - C++17 compatible compiler
%   - Eigen3 library installed
%   - OpenMP (optional, for parallelization)
%
% Installation of Eigen3:
%   macOS:   brew install eigen
%   Ubuntu:  sudo apt install libeigen3-dev
%   Windows: vcpkg install eigen3
%
% Author: Rui Yang, 2024

function compile_mex(build_type)
    if nargin < 1
        build_type = 'release';
    end

    fprintf('========================================\n');
    fprintf('AzRM MEX Compilation\n');
    fprintf('========================================\n\n');

    % Paths
    script_dir = fileparts(mfilename('fullpath'));
    cpp_dir = fullfile(script_dir, 'cpp');
    src_dir = fullfile(cpp_dir, 'src');
    inc_dir = fullfile(cpp_dir, 'include');

    % Source files
    src_files = {
        fullfile(src_dir, 'azrm_core.cpp')
        fullfile(src_dir, 'mex_interface.cpp')
    };

    % Check source files exist
    for i = 1:length(src_files)
        if ~exist(src_files{i}, 'file')
            error('Source file not found: %s', src_files{i});
        end
    end
    fprintf('Source files found.\n');

    % Find Eigen
    eigen_path = find_eigen();
    if isempty(eigen_path)
        error('Eigen3 not found. Please install Eigen3 or set EIGEN3_INCLUDE_DIR environment variable.');
    end
    fprintf('Eigen3 found: %s\n', eigen_path);

    % Compiler flags
    if ispc
        % Windows (MSVC)
        cxx_flags = 'CXXFLAGS="$CXXFLAGS /std:c++17 /O2 /openmp"';
        ld_flags = '';
    elseif ismac
        % macOS
        if strcmpi(build_type, 'debug')
            cxx_flags = 'CXXFLAGS="$CXXFLAGS -std=c++17 -O0 -g -Wall"';
        else
            cxx_flags = 'CXXFLAGS="$CXXFLAGS -std=c++17 -O3 -march=native -Wall"';
        end

        % OpenMP on macOS (via Homebrew)
        [omp_inc, omp_lib] = find_openmp_mac();
        if ~isempty(omp_inc)
            cxx_flags = [cxx_flags ' -Xpreprocessor -fopenmp -I' omp_inc];
            ld_flags = ['LDFLAGS="$LDFLAGS -L' omp_lib ' -lomp"'];
        else
            fprintf('OpenMP not found on macOS. Building without parallelization.\n');
            fprintf('Install with: brew install libomp\n');
            ld_flags = '';
        end
    else
        % Linux
        if strcmpi(build_type, 'debug')
            cxx_flags = 'CXXFLAGS="$CXXFLAGS -std=c++17 -O0 -g -Wall -fopenmp"';
        else
            cxx_flags = 'CXXFLAGS="$CXXFLAGS -std=c++17 -O3 -march=native -Wall -fopenmp"';
        end
        ld_flags = 'LDFLAGS="$LDFLAGS -fopenmp"';
    end

    % Output file
    output_name = 'azrm_mex';
    output_file = fullfile(script_dir, output_name);

    % Build mex command
    mex_args = {
        '-v'                        % Verbose
        ['-I' inc_dir]              % Include directory
        ['-I' eigen_path]           % Eigen include
        cxx_flags                   % Compiler flags
        '-output', output_file      % Output file
    };

    if ~isempty(ld_flags)
        mex_args{end+1} = ld_flags;
    end

    % Add source files
    for i = 1:length(src_files)
        mex_args{end+1} = src_files{i};
    end

    % Print command
    fprintf('\nMEX command:\n');
    fprintf('  mex %s\n\n', strjoin(mex_args, ' '));

    % Compile
    fprintf('Compiling...\n');
    try
        mex(mex_args{:});
        fprintf('\n========================================\n');
        fprintf('Compilation successful!\n');
        fprintf('MEX file: %s.%s\n', output_file, mexext);
        fprintf('========================================\n\n');

        % Quick test
        fprintf('Running quick test...\n');
        quick_test(output_file);

    catch ME
        fprintf('\n========================================\n');
        fprintf('Compilation failed!\n');
        fprintf('Error: %s\n', ME.message);
        fprintf('========================================\n\n');
        fprintf('Troubleshooting:\n');
        fprintf('1. Make sure you have a C++17 compatible compiler\n');
        fprintf('2. Run "mex -setup C++" to configure compiler\n');
        fprintf('3. Install Eigen3: brew install eigen (macOS)\n');
        fprintf('4. For OpenMP on macOS: brew install libomp\n');
        rethrow(ME);
    end
end

function eigen_path = find_eigen()
    % Check environment variable
    eigen_path = getenv('EIGEN3_INCLUDE_DIR');
    if ~isempty(eigen_path) && exist(fullfile(eigen_path, 'Eigen', 'Core'), 'file')
        return;
    end

    % Common paths
    if ismac
        paths = {
            '/opt/homebrew/include/eigen3'
            '/usr/local/include/eigen3'
            '/opt/local/include/eigen3'
        };
    elseif isunix
        paths = {
            '/usr/include/eigen3'
            '/usr/local/include/eigen3'
        };
    else  % Windows
        paths = {
            'C:\vcpkg\installed\x64-windows\include\eigen3'
            'C:\Eigen3\include'
        };
    end

    for i = 1:length(paths)
        if exist(fullfile(paths{i}, 'Eigen', 'Core'), 'file')
            eigen_path = paths{i};
            return;
        end
    end

    eigen_path = '';
end

function [omp_inc, omp_lib] = find_openmp_mac()
    omp_inc = '';
    omp_lib = '';

    % Homebrew paths
    paths = {
        '/opt/homebrew/opt/libomp'
        '/usr/local/opt/libomp'
    };

    for i = 1:length(paths)
        if exist(fullfile(paths{i}, 'include', 'omp.h'), 'file')
            omp_inc = fullfile(paths{i}, 'include');
            omp_lib = fullfile(paths{i}, 'lib');
            return;
        end
    end
end

function quick_test(mex_file)
    % Add to path
    addpath(fileparts(mex_file));

    % Simple two-layer model
    nL = 2;
    thickness = [0.1; 0.5];
    den = [2.2; 2.5];

    % VTI stiffness
    c_layer = zeros(6, 6, nL);
    c_layer(:,:,1) = Stiffness_matrix_VTI(2.5, 1.2, 2.2, 0.05, 0.02, 0.03);
    c_layer(:,:,2) = Stiffness_matrix_VTI(3.5, 1.8, 2.5, 0.08, 0.04, 0.05);

    freq = 30;
    angles = [0, 10, 20, 30];
    phi = 0;

    % Call MEX
    [Rpp, ~, ~] = azrm_mex(freq, angles, thickness, den, c_layer, phi);

    fprintf('Quick test results:\n');
    fprintf('  Rpp at 0 deg:  %.6f + %.6fi\n', real(Rpp(1)), imag(Rpp(1)));
    fprintf('  Rpp at 30 deg: %.6f + %.6fi\n', real(Rpp(4)), imag(Rpp(4)));

    % Validate
    if all(isfinite(Rpp)) && all(abs(Rpp) < 1)
        fprintf('Quick test: PASSED\n');
    else
        fprintf('Quick test: FAILED\n');
    end
end
