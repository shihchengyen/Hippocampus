function filename = hash_filename(Args)
    hash_input = '';
    for i = 1:length(Args.DataCheckArgs)
        key = Args.DataCheckArgs{i};
        val = Args.(key);
        if isnumeric(val)
            valStr = mat2str(val);
        elseif ischar(val)
            valStr = val;
        else
            error(['Unsupported argument type: ' key]);
        end
        hash_input = [hash_input '_'  valStr];
    end
    hashString = javaHash(hash_input, 'SHA-256');
    filename = sprintf('%s_%s.mat', Args.classname, hashString(1:16));  % Shorten to 16 chars

function hash = javaHash(inputStr, algorithm)
    % Convert string to bytes (uint8)
    data = uint8(inputStr);
    % Create MessageDigest object for the specified algorithm (e.g., SHA-256)
    md = java.security.MessageDigest.getInstance(algorithm);
    % Update the MessageDigest with the data
    md.update(data);
    % Get the resulting hash (as a byte array)
    hashBytes = md.digest();
    % Convert hash to a hex string
    hash = '';
    for i = 1:length(hashBytes)
        hash = [hash, sprintf('%02x', hashBytes(i))];  % Format each byte as 2-digit hex
    end

