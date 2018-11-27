classdef Flocking_App_final < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        UIAxes                          matlab.ui.control.UIAxes
        KEditFieldLabel                 matlab.ui.control.Label
        KEditField                      matlab.ui.control.NumericEditField
        SigmaEditFieldLabel             matlab.ui.control.Label
        SigmaEditField                  matlab.ui.control.NumericEditField
        BetaEditFieldLabel              matlab.ui.control.Label
        BetaEditField                   matlab.ui.control.NumericEditField
        ofAgentsEditFieldLabel          matlab.ui.control.Label
        ofAgentsEditField               matlab.ui.control.NumericEditField
        LeaderTrajectoryLabel           matlab.ui.control.Label
        UITable                         matlab.ui.control.Table
        PlotButton                      matlab.ui.control.Button
        PlayButton                      matlab.ui.control.Button
        StopButton                      matlab.ui.control.Button
        ResetButton                     matlab.ui.control.Button
        XEditFieldLabel                 matlab.ui.control.Label
        XEditField                      matlab.ui.control.EditField
        YEditFieldLabel                 matlab.ui.control.Label
        YEditField                      matlab.ui.control.EditField
        MaxIterationsEditFieldLabel     matlab.ui.control.Label
        MaxIterationsEditField          matlab.ui.control.NumericEditField
        CurrentIterationEditFieldLabel  matlab.ui.control.Label
        CurrentIterationEditField       matlab.ui.control.NumericEditField
        LeaderLabel                     matlab.ui.control.Label
    end

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
    %%%% Fill table with random values during startup %%%%
    
    %If you want to specifically define agent starting positions, alter the c matrix here. 
    %Similar code should also be altered in the reset callback, 
    %so pressing reset brings your agents back to the correct positions
    
    n = app.ofAgentsEditField.Value;
    for i = 1:n
    c(i,1) = 0;
    c(i,2) = 0;
    c(i,3) = 0;
    c(i,4) = 0;
    c(i,5) = 0;
    end
    
    app.UITable.Data = c;
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
%This is the main code that executes the flocking algorithm when the play button is pressed
            
app.CurrentIterationEditField.Value = 1;
max_iterations = app.MaxIterationsEditField.Value;

%Get the number of agents
n = app.ofAgentsEditField.Value;

%Get the initial locations and velocities of the agents and store them in
%Hx and Hv respectively
data=app.UITable.Data;
Hx=data(:,2:3);
Hv=data(:,4:5);

%Get the K, sigma and beta values for use in the function for defining the
%adjacency matrix
K = app.KEditField.Value;
sigma = app.SigmaEditField.Value;
beta = app.BetaEditField.Value;

%Initially set the size of the adjacency matrix
A=zeros(n);

%Choose time step. YOU CAN CHANGE THIS. 
deltat=0.05; 

time=1;

arena_dimensions = 5;

%set up symbolic variable used when interpreting trajectories
syms t 

%Get initial position of leader and x,y trajectories of leader
x_pos(1) = data(1,1);
y_pos(1) = data(1,2);
func_x = str2sym(app.XEditField.Value);
func_y = str2sym(app.YEditField.Value);

%Evaluate the positions of the leader over its trajectory.
%The leader will start at its initial position in the data table
%and will execute one time step with its initial velocity before
%following the trajectory from its current position. 
%(The leader is agent 1 in the data table)

if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '') %IF we are given a leader trajectory, use trajectory for leader
    for k = 2:app.MaxIterationsEditField.Value
        t = k;
            x_pos(k) = subs(func_x) + x_pos(1);
            y_pos(k) = subs(func_y) + y_pos(1);
    end 
end 

%--------------------------------------------------------------------------
%Input your trigger sequence here to determine on which time iterations your
%agents will update their velocities
%Call your sequence trigger and make it a vector

trigger=zeros(1,app.MaxIterationsEditField.Value); %%%% Default Trigger (always fires) %%%%
for i=1:app.MaxIterationsEditField.Value
    if mod(i,5) == 0
        trigger(i) = 1;
    end
end
%--------------------------------------------------------------------------

%Execution Loop 

while ((app.CurrentIterationEditField.Value < (max_iterations)-1) && (app.CurrentIterationEditField.Value > 0))
    %Define the entries of A by the distances of the nodes
    for i=1:n
        for j=1:n
            dist=pdist([Hx(i,1) Hx(i,2);Hx(j,1) Hx(j,2)]);
            A(i,j)=K/((sigma^2+dist^2)^beta);
        end
    end
    
if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
    %This executes if we are given a leader. 
%   Store the x and y positions of the leader in a single vector (pos)
    for i=1:size(x_pos,2)-1
        pos(2*i-1)=x_pos(i);
        pos(2*i)=y_pos(i);
    end
    
    %Calculate the velocity of leader
    for i=1:((time+1)/2)+1
        x_vel(i)=x_pos(i+1)-x_pos(i);
        y_vel(i)=y_pos(i+1)-y_pos(i);
    end
    
    %Store the x and y velocity in a single vector (vel)
    for i=1:size(x_vel,2)
        vel(2*i-1)=x_vel(i);
        vel(2*i)=y_vel(i);
    end 
end 

    %Update positions using the algorithm
    Hx(:,time+2)=Hx(:,time)+deltat*Hv(:,time);
    Hx(:,time+3)=Hx(:,time+1)+deltat*Hv(:,time+1);
    
    %Compute Laplacian Matrix L
    w=A*ones(n,1);
    D=diag(w);
    L=D-A;
    if Hv(1,time+1) > Hv(2,time+1)
        L(2,1) = 1/L(2,1);
    end
    if Hv(2,time+1) > Hv(3,time+1)
        L(2,3) = 1/L(2,3);
    end
    if Hv(1,time+1) > Hv(3,time+1)
        L(3,1) = 1/L(3,1);
    end
    if Hv(2,time+1) > Hv(3,time+1)
        L(3,2) = 1/L(3,2);
    end
    if trigger((time+1)/2)==1
        Hv(:,time+2)=Hv(:,time)-deltat.*L*Hv(:,time);
        Hv(:,time+3)=Hv(:,time+1)-deltat.*L*Hv(:,time+1);
    else
        Hv(:,time+2)=Hv(:,time);
        Hv(:,time+3)=Hv(:,time+1);
    end
    
    %Reset the leader's (ie first node) velocity to what was
    %calculated above
    if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
     Hv(1,:)=vel;
    end 
    
    %Increment the time by 2 since the Hx=[x_position,y_position] and
    % Hv=[x_vel,y_vel]
    app.CurrentIterationEditField.Value = app.CurrentIterationEditField.Value + 1;
    time=time+2;
    
    %Set the current table data in the GUI
    tabledata=[Hx(:,[time,time+1]) Hv(:,[time,time+1])];
    app.UITable.Data(:, 2:5) = tabledata;

    %----------------------------------------------------------------------
    %Plot the nodes
    %----------------------------------------------------------------------
    app.UIAxes.ColorOrderIndex = 1;
    if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
        plot(app.UIAxes,Hx(1,time),Hx(1,time+1),'Marker','o','MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') %want to fill the leader in with a blue face
    else 
        plot(app.UIAxes,Hx(1,time),Hx(1,time+1),'Marker','o','MarkerSize',6)
    end 
    
    hold(app.UIAxes, 'on')
    for i=2:n
        plot(app.UIAxes,Hx(i,time),Hx(i,time+1),'Marker','o','MarkerSize',6)
        hold(app.UIAxes, 'on')
    end
    if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
     text(app.UIAxes, Hx(1,time),Hx(1,time+1), '\leftarrow Leader', 'color', 'k', 'fontweight', 'bold')
    end 
    %Construct x and y arrays for use in plotting the agent trajectories
    for i=1:n
        for j=1:2:time
            x(i,(j+1)/2)=Hx(i,j);
        end
        for s=2:2:time+1
            y(i,s/2)=Hx(i,s);
        end
    end
    
    app.UIAxes.ColorOrderIndex = 1;
    %Plot the agent trajectories
    for i=1:n
        plot(app.UIAxes,x(i,:),y(i,:));
        hold(app.UIAxes, 'on')
        
            
        %Plot agent numbers
        text(app.UIAxes, x(i,app.CurrentIterationEditField.Value), y(i,app.CurrentIterationEditField.Value), sprintf('%d', i), 'color', 'k', ...
                        'fontweight', 'bold');
    end
    
    hold(app.UIAxes, 'off')
    pause(0.05)
    
end 

%Write history to excel spreadsheet called agentData. Column 1
%contains agent 1's x positions, column 2 contains agent 2's y
%positions, column 3 contains agent 2's x positions, and so on. The first
%row contains initial positions, the following rows containing the next
%position, and so on, until the final position.

Ax = zeros(n,size(Hx,2)/2);
Ay = zeros(n,size(Hx,2)/2);
j = 1;

for i = 1 : 2 : size(Hx,2)-1
    Ax(:,j) = Hx(:,i);
    Ay(:,j) = Hx(:,i+1);
    j = j+1;
end

Ax = transpose(Ax);
Ay = transpose(Ay); 
H = horzcat(Ax(:,1),Ay(:,1));

for i = 2 : n
    H = horzcat(H,Ax(:,i),Ay(:,i));
end

%------------------- EXCEL SHEETS -------------------% 

min_value = min(min(H));
disp(min_value);

if min_value < 0
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            H1(i,j) = H(i,j) - min_value; % Get everything to be positive
        end
    end
else
H1 = H;
end

max_value = max(max(H1));
disp(max_value);
H1 = H1*(arena_dimensions/max_value); %Resizing
delete agentData.xlsx
xlswrite('agentData.xlsx', H, 1);
xlswrite('agentData.xlsx', H1, 2);
save('agentData.mat','H');
    
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
    %Stop plotting when the stop button is pushed
    %Achieved by bumping the current interation up to the final iteration
        
    app.CurrentIterationEditField.Value = app.MaxIterationsEditField.Value + 1;
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
   %Clear figure and reset data when reset button is pressed
        
    app.CurrentIterationEditField.Value = 0;
    cla(app.UIAxes);
    
    for i = 1:app.ofAgentsEditField.Value
    c(i,1) = i;
    c(i,2) = randi([-10 10]);
    c(i,3) = randi([-10 10]);
    c(i,4) = randi([-10 10]);
    c(i,5) = randi([-10 10]);
    end
    
    app.UITable.Data = c;
        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
    %Plot table data when the plot button is pushed
      data = app.UITable.Data(:,2:3);  
      app.UIAxes.ColorOrderIndex = 1;
      
    if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
        plot(app.UIAxes,data(1,1),data(1,2),'Marker','o','MarkerSize',6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
    else 
        plot(app.UIAxes,data(1,1),data(1,2),'Marker','o','MarkerSize',6)
    end 
    n = app.ofAgentsEditField.Value;
    hold(app.UIAxes, 'on')
    for i=2:n 
    plot(app.UIAxes,data(i,1),data(i,2),'Marker','o','MarkerSize',6)
    hold(app.UIAxes, 'on')
    
     %Plot agent numbers
         text(app.UIAxes, data(i,1), data(i,2), sprintf('%d', i), 'color', 'k', ...
                         'fontweight', 'bold');
    end
    
    if ~strcmpi(app.XEditField.Value,'') && ~strcmpi(app.YEditField.Value, '')
        text(app.UIAxes, data(1,1),data(1,2), '\leftarrow Leader', 'color', 'k', 'fontweight', 'bold')
    else 
        text(app.UIAxes, data(1,1), data(1,2), '1', 'color', 'k', 'fontweight', 'bold')
    end 
    
    hold(app.UIAxes, 'off')
        end

        % Value changed function: ofAgentsEditField
        function ofAgentsEditFieldValueChanged(app, event)
    %Update table size when the number of agents is changed. 
    
    n = app.ofAgentsEditField.Value;
    
       for i = 1:n
        c(i,1) = i;
        c(i,2) = randi([-10 10]);
        c(i,3) = randi([-10 10]);
        c(i,4) = randi([-10 10]);
        c(i,5) = randi([-10 10]);
       end
    
    app.UITable.Data = c;
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 753 565];
            app.UIFigure.Name = 'UI Figure';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Agent Position')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.Position = [8 18 400 541];

            % Create KEditFieldLabel
            app.KEditFieldLabel = uilabel(app.UIFigure);
            app.KEditFieldLabel.HorizontalAlignment = 'right';
            app.KEditFieldLabel.Position = [622 523 25 22];
            app.KEditFieldLabel.Text = 'K';

            % Create KEditField
            app.KEditField = uieditfield(app.UIFigure, 'numeric');
            app.KEditField.Position = [657 523 45 22];
            app.KEditField.Value = 1;

            % Create SigmaEditFieldLabel
            app.SigmaEditFieldLabel = uilabel(app.UIFigure);
            app.SigmaEditFieldLabel.HorizontalAlignment = 'right';
            app.SigmaEditFieldLabel.Position = [603 489 41 22];
            app.SigmaEditFieldLabel.Text = 'Sigma';

            % Create SigmaEditField
            app.SigmaEditField = uieditfield(app.UIFigure, 'numeric');
            app.SigmaEditField.Position = [657 489 44 22];
            app.SigmaEditField.Value = 1;

            % Create BetaEditFieldLabel
            app.BetaEditFieldLabel = uilabel(app.UIFigure);
            app.BetaEditFieldLabel.HorizontalAlignment = 'right';
            app.BetaEditFieldLabel.Position = [614 456 30 22];
            app.BetaEditFieldLabel.Text = 'Beta';

            % Create BetaEditField
            app.BetaEditField = uieditfield(app.UIFigure, 'numeric');
            app.BetaEditField.Position = [657 456 44 22];
            app.BetaEditField.Value = 0.25;

            % Create ofAgentsEditFieldLabel
            app.ofAgentsEditFieldLabel = uilabel(app.UIFigure);
            app.ofAgentsEditFieldLabel.HorizontalAlignment = 'right';
            app.ofAgentsEditFieldLabel.Position = [578 421 66 22];
            app.ofAgentsEditFieldLabel.Text = '# of Agents';

            % Create ofAgentsEditField
            app.ofAgentsEditField = uieditfield(app.UIFigure, 'numeric');
            app.ofAgentsEditField.ValueChangedFcn = createCallbackFcn(app, @ofAgentsEditFieldValueChanged, true);
            app.ofAgentsEditField.Position = [657 421 44 22];
            app.ofAgentsEditField.Value = 10;

            % Create LeaderTrajectoryLabel
            app.LeaderTrajectoryLabel = uilabel(app.UIFigure);
            app.LeaderTrajectoryLabel.Position = [447 523 100 22];
            app.LeaderTrajectoryLabel.Text = 'Leader Trajectory';

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'Agent'; 'X'; 'Y'; 'Vx'; 'Vy'};
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = [false true true true true];
            app.UITable.Position = [479 154 252 247];

            % Create PlotButton
            app.PlotButton = uibutton(app.UIFigure, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Position = [479 108 100 24];
            app.PlotButton.Text = 'Plot';

            % Create PlayButton
            app.PlayButton = uibutton(app.UIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.Position = [602 108 100 24];
            app.PlayButton.Text = 'Play';

            % Create StopButton
            app.StopButton = uibutton(app.UIFigure, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Position = [479 65 100 24];
            app.StopButton.Text = 'Stop';

            % Create ResetButton
            app.ResetButton = uibutton(app.UIFigure, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [601 65 100 24];
            app.ResetButton.Text = 'Reset';

            % Create XEditFieldLabel
            app.XEditFieldLabel = uilabel(app.UIFigure);
            app.XEditFieldLabel.HorizontalAlignment = 'right';
            app.XEditFieldLabel.Position = [407 489 25 22];
            app.XEditFieldLabel.Text = 'X';

            % Create XEditField
            app.XEditField = uieditfield(app.UIFigure, 'text');
            app.XEditField.Position = [447 489 100 22];

            % Create YEditFieldLabel
            app.YEditFieldLabel = uilabel(app.UIFigure);
            app.YEditFieldLabel.HorizontalAlignment = 'right';
            app.YEditFieldLabel.Position = [407 456 25 22];
            app.YEditFieldLabel.Text = 'Y';

            % Create YEditField
            app.YEditField = uieditfield(app.UIFigure, 'text');
            app.YEditField.Position = [447 456 100 22];

            % Create MaxIterationsEditFieldLabel
            app.MaxIterationsEditFieldLabel = uilabel(app.UIFigure);
            app.MaxIterationsEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxIterationsEditFieldLabel.Position = [410 421 80 22];
            app.MaxIterationsEditFieldLabel.Text = 'Max Iterations';

            % Create MaxIterationsEditField
            app.MaxIterationsEditField = uieditfield(app.UIFigure, 'numeric');
            app.MaxIterationsEditField.Position = [500 421 47 22];
            app.MaxIterationsEditField.Value = 50;

            % Create CurrentIterationEditFieldLabel
            app.CurrentIterationEditFieldLabel = uilabel(app.UIFigure);
            app.CurrentIterationEditFieldLabel.HorizontalAlignment = 'right';
            app.CurrentIterationEditFieldLabel.Position = [493 18 93 22];
            app.CurrentIterationEditFieldLabel.Text = 'Current Iteration';

            % Create CurrentIterationEditField
            app.CurrentIterationEditField = uieditfield(app.UIFigure, 'numeric');
            app.CurrentIterationEditField.Position = [601 18 100 22];

            % Create LeaderLabel
            app.LeaderLabel = uilabel(app.UIFigure);
            app.LeaderLabel.Position = [424 359 52 22];
            app.LeaderLabel.Text = '(Leader)';
        end
    end

    methods (Access = public)

        % Construct app
        function app = Flocking_App_final

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end