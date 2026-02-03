```python
import json
import hashlib
import os
from typing import Dict, Any

# Define the path for the users data file
_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
_USERS_FILE_PATH = os.path.join(_DATA_DIR, 'users.json')

def _ensure_data_dir_exists():
    """Ensures the data directory exists."""
    os.makedirs(_DATA_DIR, exist_ok=True)

def _load_users() -> Dict[str, Any]:
    """
    Loads user data from the JSON file.

    Returns:
        Dict[str, Any]: A dictionary where keys are lowercase usernames
                        and values are user data dictionaries.
    """
    _ensure_data_dir_exists()
    if not os.path.exists(_USERS_FILE_PATH):
        return {}
    try:
        with open(_USERS_FILE_PATH, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError:
        print(f"Warning: Corrupted or empty users.json at {_USERS_FILE_PATH}. Starting with an empty user database.")
        return {}
    except Exception as e:
        print(f"Error loading users from {_USERS_FILE_PATH}: {e}")
        return {}

def _save_users(users: Dict[str, Any]):
    """
    Saves user data to the JSON file.

    Args:
        users (Dict[str, Any]): The dictionary of user data to save.
    """
    _ensure_data_dir_exists()
    try:
        with open(_USERS_FILE_PATH, 'w') as f:
            json.dump(users, f, indent=4)
    except Exception as e:
        print(f"Error saving users to {_USERS_FILE_PATH}: {e}")

def _generate_salt() -> str:
    """
    Generates a random salt for password hashing.

    Returns:
        str: A hexadecimal string representing the salt.
    """
    return os.urandom(16).hex()

def _hash_password(password: str, salt: str) -> str:
    """
    Hashes a password using SHA256 and a given salt.

    Args:
        password (str): The plain-text password.
        salt (str): The salt to use for hashing.

    Returns:
        str: The SHA256 hexadecimal digest of the salted password.
    """
    salted_password = (password + salt).encode('utf-8')
    return hashlib.sha256(salted_password).hexdigest()

def register_user(username: str, password: str) -> bool:
    """
    Registers a new user in the system.

    Args:
        username (str): The desired username (case-insensitive for uniqueness check).
        password (str): The desired password.

    Returns:
        bool: True if registration was successful, False if the username already exists
              or if the username/password is invalid (e.g., empty).
    """
    if not username or not password:
        print("Registration failed: Username and password cannot be empty.")
        return False

    users = _load_users()
    
    # Use lowercase username for internal lookup to ensure case-insensitive uniqueness
    normalized_username = username.lower()

    if normalized_username in users:
        print(f"Registration failed: Username '{username}' already exists.")
        return False

    salt = _generate_salt()
    hashed_password = _hash_password(password, salt)

    users[normalized_username] = {
        "username": username, # Store original casing for display purposes if needed
        "password_hash": hashed_password,
        "salt": salt
    }
    _save_users(users)
    print(f"User '{username}' registered successfully.")
    return True

def login_user(username: str, password: str) -> bool:
    """
    Authenticates a user against the stored credentials.

    Args:
        username (str): The username provided for login (case-insensitive).
        password (str): The password provided for login.

    Returns:
        bool: True if login was successful, False otherwise (username not found or incorrect password).
    """
    if not username or not password:
        return False

    users = _load_users()
    
    # Use lowercase username for internal lookup
    normalized_username = username.lower()
    user_data = users.get(normalized_username)

    if not user_data:
        # For security, avoid distinguishing between "username not found" and "incorrect password"
        return False

    stored_hash = user_data["password_hash"]
    stored_salt = user_data["salt"]

    provided_password_hash = _hash_password(password, stored_salt)

    if provided_password_hash == stored_hash:
        print(f"Login successful for '{user_data['username']}'.")
        return True
    else:
        # For security, avoid distinguishing between "username not found" and "incorrect password"
        return False

if __name__ == "__main__":
    print("--- User Authentication Module Test ---")

    # Clean up previous data for a fresh run
    if os.path.exists(_USERS_FILE_PATH):
        os.remove(_USERS_FILE_PATH)
        print(f"Cleaned up existing '{_USERS_FILE_PATH}' for test run.")

    print("\n--- Registering Users ---")
    register_user("testuser", "password123")  # Success
    register_user("TestUser", "password123")  # Should fail (case-insensitive duplicate)
    register_user("admin", "adminpass")       # Success
    register_user("", "empty")                # Should fail (empty username)
    register_user("valid", "")                # Should fail (empty password)

    print("\n--- Logging In Users ---")
    print(f"Login 'testuser' with 'password123': {login_user('testuser', 'password123')}") # Success
    print(f"Login 'TestUser' with 'password123': {login_user('TestUser', 'password123')}") # Success (case-insensitive)
    print(f"Login 'testuser' with 'wrongpass': {login_user('testuser', 'wrongpass')}")     # Failure
    print(f"Login 'nonexistent' with 'anypass': {login_user('nonexistent', 'anypass')}")   # Failure
    print(f"Login 'admin' with 'adminpass': {login_user('admin', 'adminpass')}")           # Success
    print(f"Login 'admin' with 'admin_pass': {login_user('admin', 'admin_pass')}")         # Failure (wrong password)

    print(f"\n--- Current Content of '{_USERS_FILE_PATH}' ---")
    try:
        with open(_USERS_FILE_PATH, 'r') as f:
            print(f.read())
    except FileNotFoundError:
        print("Users file does not exist (unexpected after registrations).")
    except Exception as e:
        print(f"Error reading users file: {e}")

```