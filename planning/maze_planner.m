
function [planner, planner2, start, goal, waypoints,trajA,trajAS,trajR,trajRS]= maze_planner(filename, map, inflate)
%MAZE_PLANNER - Path planning for maze navigation using Hybrid A* and RRT*
%   Computes a collision-free reference path from START to GOAL positions
%   defined in an XML world file using Hybrid A* and RRT* algorithm on an occupancy
%   map.
%
% Input Arguments:
%   filename (string)    - Path to XML world file containing START and GOAL markers
%   map (H x W x C)      - Binary or grayscale image map (only first channel used)
%   inflate (double)     - Obstacle inflation radius [m] (default: 0.3)
%
% Output Arguments:
%   planner (plannerHybridAStar) -  Configured Hybrid A* planner object
%   planner2 (plannerRRT*) -        Configured RRT* planner object
%   start (1x3)                  - Start pose [x, y, theta] in world frame [m, m, rad]
%   goal (1x3)                   - Goal pose [x, y, theta] in world frame [m, m, rad]
%   waypoints (Nx3)              - Waypoints along the planned path [x, y, theta]
%   trajA (navPath)             - Planned path object containing states and
%   directions from Planner Hybrid A*  
%   trajAS (navPath)            - Smoth Planned path object containing states and
%   directions from Planner Hybrid A* 
%   trajR (navPath)             - Direct Planned path object containing states and
%   directions from Planner RRT*
%   trajRS (navPath)            - Direct Smoth Planned path object containing
%   states and directions from Planner RRT*

% Example:
%   map = imread('maze.png');
%   [planner, start, goal, path,trajA,trajAS,trajR,trajRS] = maze_planner('world.xml', map);
%   show(path);

arguments
    filename {mustBeTextScalar}
    map (:,:,:) {mustBeNumeric}
    inflate (1,1) double {mustBeNonnegative, mustBeFinite} = 0.3
end

%% Parse World XML File
worldxml = readstruct(filename, "FileType", "xml");

waypoints = [2 4 0 0];  % add start

it = 2;

% Waypoints based on statues
    for i = 1:size(worldxml.World.MeshPart, 2)
        part = worldxml.World.MeshPart(i);
        if isa(part.position, 'string')
            coords = double(split(extractBetween(part.position, "{", "}"), ","));
        end
        if isa(part.orientation, 'string')
            orient = split(extractBetween(part.orientation, "{", "}"), ",");
        else
            orient = [0, 0, 0];
        end
        if part.nameAttribute.startsWith("death")
            LM = part.nameAttribute.split('-');LM = LM(2);
            waypointNum2 = LM.replace("LM","");
            ang = rad2deg(double(orient(3)));
            if ang < 0
                ang = 360 + ang;
            end
            if ang < 15 || ang > 345 % Looking up
                coords(2) = coords(2)+2; 
            elseif ang < 195 && ang > 165 % Looking down
                coords(2) = coords(2)-2;
            elseif ang > 255 && ang < 285 % Looking right
                coords(1) = coords(1)+2; 
            elseif ang > 75 && ang < 105 % Looking left
                coords(1) = coords(1)-2;
            end
            waypoints(it,:) = [coords(1),coords(2),0,double(waypointNum2)];
            it=it+1;
        end
    end

waypoints=[waypoints;[30 56 0 0]];   %% add goal 
waypoints(:,4) = [];        % clear it waypointNum orient coords;

%% Configure State Space and Validator
% SE(2) state space: [x, y, theta]
state_space = stateSpaceSE2;
state_validator = validatorOccupancyMap(state_space);

% Create binary occupancy map from image (rotate 180° for coordinate alignment)
occupancy_map = binaryOccupancyMap(imrotate(map(:,:,1), 180), "Resolution", 10);
occupancy_map.inflate(inflate);

% Configure validator
state_validator.Map = occupancy_map;
state_validator.ValidationDistance = 0.1;

% Set state bounds based on map limits
state_space.StateBounds = [
    occupancy_map.XWorldLimits;
    occupancy_map.YWorldLimits;
    [-pi pi]
];

%% Compute Reference Trajectory

% Hybrid A* planner for car-like robots

planner = plannerHybridAStar(state_validator,"MinTurningRadius",0.65,"DirectionSwitchingCost",10);

start=[2 4 0];
goal=[30 56 0];

[x, ~]=size(waypoints);
trajA=[];
trajAS=[];                % Smooth path

for i=1:x-1
    refpath = plan(planner,waypoints(i,:),waypoints(i+1,:));
    trajAux = refpath.States;

    trajAux(end,:)=[];
    trajA=[trajA;trajAux];

end

trajAS(:,1) = smoothdata(trajA(:,1), 'gaussian', 4);
trajAS(:,2) = smoothdata(trajA(:,2), 'gaussian', 4);
trajAS(:,3) = trajA(:,3);

% RRT* Planner

planner2 = plannerRRTStar(state_space,state_validator,"MaxConnectionDistance",2.5,"BallRadiusConstant",100);

refpath2 = plan(planner2,start,goal);

trajR=refpath2.States;

trajRS(:,1) = smoothdata(trajR(:,1), 'gaussian', 4);
trajRS(:,2) = smoothdata(trajR(:,2), 'gaussian', 4);
trajRS(:,3) = trajR(:,3);

end