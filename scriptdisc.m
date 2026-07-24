% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
%
% Reference:
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
function scriptdisc(filename)
% scriptdisc  Print a console disclaimer/citation notice for filename.
%
%   scriptdisc(filename)
%
%   Not called by InstanceSpace.build()/explore() or the buildIS/exploreIS
%   wrappers as of the Phase 7 refactor; kept for any external callers
%   still invoking it directly. Superseded by the LICENSE file and the
%   per-file SPDX headers (spec Section 11); slated for retirement in a
%   future repository-hygiene pass (spec Phase 10).

disp('-------------------------------------------------------------------------');
disp(' ');
disp([' ''' filename ''' ']);
disp(' ');
disp(' Code by: Mario Andres Mu�oz Acosta');
disp('          School of Mathematics and Statistics');
disp('          The University of Melbourne');
disp('          Australia');
disp('          2019');
disp(' ');
disp(' Copyright: Mario A. Mu�oz');
disp(' ');
disp('-------------------------------------------------------------------------');
disp(' ');
disp(' If using this software, please cite as: ');
disp(' ');
disp([' Mario Andr�s Mu�oz, & Kate Smith-Miles. ' ...
      ' andremun/InstanceSpace: February 2021 Update (Version v0.2-beta). '...
      ' Zenodo. http://doi.org/10.5281/zenodo.4521336']);
disp(' ');
disp('-------------------------------------------------------------------------');
disp(' ');
disp(' DISCLAIMER:');
disp(' ');
disp(' THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY');
disp(' APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT');
disp(' HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY');
disp(' OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,');
disp(' THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR');
disp(' PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM');
disp(' IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF');
disp(' ALL NECESSARY SERVICING, REPAIR OR CORRECTION.');
disp(' ');
disp(' IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING');
disp(' WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS');
disp(' THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY');
disp(' GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE');
disp(' USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF');
disp(' DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD');
disp(' PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),');
disp(' EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF');
disp(' SUCH DAMAGES.');
disp(' ');

end