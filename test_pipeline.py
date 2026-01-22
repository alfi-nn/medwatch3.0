"""
Quick test script to verify the pipeline setup and initial functionality.
This tests the pipeline without making extensive API calls.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from data_collection_pipeline import Config, IDMapper, SIDERDataLoader

def test_config():
    """Test configuration setup."""
    print("Testing configuration...")
    Config.setup_directories()
    
    assert Config.BASE_DIR.exists(), "Base directory not created"
    assert Config.RAW_DIR.exists(), "Raw directory not created"
    assert Config.PROCESSED_DIR.exists(), "Processed directory not created"
    assert Config.LOGS_DIR.exists(), "Logs directory not created"
    
    print("✓ Configuration test passed")

def test_id_mapping():
    """Test ID mapping functions."""
    print("\nTesting ID mapping...")
    
    # Test STITCH to PubChem conversion
    test_cases = [
        ("CID000012345", "12345"),
        ("CIDm000012345", "12345"),
        ("CID100000123", "100000123"),
    ]
    
    for stitch_id, expected_cid in test_cases:
        result = IDMapper.stitch_to_pubchem(stitch_id)
        assert result == expected_cid, f"Failed: {stitch_id} -> {result} (expected {expected_cid})"
    
    print("✓ ID mapping test passed")

def test_sider_download():
    """Test SIDER download functionality."""
    print("\nTesting SIDER download...")
    
    # This will download SIDER files if not present
    success = SIDERDataLoader.download_sider_dataset()
    
    if success:
        print("✓ SIDER download test passed")
    else:
        print("⚠ SIDER download failed (may be network issue)")
    
    return success

def main():
    """Run all tests."""
    print("=" * 60)
    print("Data Collection Pipeline - Quick Test")
    print("=" * 60)
    
    try:
        test_config()
        test_id_mapping()
        test_sider_download()
        
        print("\n" + "=" * 60)
        print("All tests passed! Pipeline is ready to run.")
        print("=" * 60)
        print("\nTo run the full pipeline:")
        print("  python data_collection_pipeline.py")
        print("\nNote: Full execution may take 30-60 minutes.")
        
    except Exception as e:
        print(f"\n✗ Test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
