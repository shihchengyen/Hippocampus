function shortestPathLength = bfs(adjList, startNode, targetNode)
    visited = false(size(adjList));
    queue = [startNode];
    distances = zeros(size(adjList));
    distances(:) = Inf;
    distances(startNode) = 0;
    
    while ~isempty(queue)
        currentNode = queue(1);
        queue(1) = [];
        
        for neighbor = adjList{currentNode}
            if ~visited(neighbor)
                visited(neighbor) = true;
                queue = [queue, neighbor];
                distances(neighbor) = distances(currentNode) + 1;
            end
        end
    end
    
    shortestPathLength = distances(targetNode);
end
