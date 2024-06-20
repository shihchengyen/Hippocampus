function shortestPathLength = bfs(graph, startNode, targetNode)
    if iscell(graph)
        shortestPathLength = bfsAdjLst(graph, startNode, targetNode);
    elseif isnumeric(graph)
        shortestPathLength = bfsAdjMat(graph, startNode, targetNode);
    else
        shortestPathLength = NaN;
    end
end

function shortestPathLength = bfsAdjLst(graph, source, target)
    visited = false(size(graph));
    distances = nan(size(graph));
    distances(source) = 0;
    queue = [source]; 
    
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];
        
        if current == target
            break;
        end
        
        for neighbour = graph{current}
            if ~visited(neighbour)
                visited(neighbour) = true;
                distances(neighbour) = distances(current) + 1;
                queue = [queue, neighbour];
            end
        end
    end
    
    shortestPathLength = distances(target);
end


function shortestPathLength = bfsAdjMat(graph, source, target)
    visited = false(size(graph));
    distances = nan(size(graph));
    distances(source) = 0;
    queue = [source]; 
    
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];
        
        if current == target
            break;
        end
        
        neighbours = find(graph(current, :));
        for neighbour = neighbours
            if ~visited(neighbour)
                visited(neighbour) = true;
                distances(neighbour) = distances(current) + 1;
                queue = [queue, neighbour];
            end
        end
    end
    
    shortestPathLength = distances(target);
end
