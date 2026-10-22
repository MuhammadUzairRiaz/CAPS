// Release manifests and the update decision (caps/update.hpp): parsing caps-release/1, versions first (0.1.0 → 0.2.0
// is an update, a lower release is not), then with equal versions the commit and its date, the opt-out variable, and
// the download URL for a platform. No network: the feed is a local file.
#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "caps/update.hpp"

using namespace caps;

namespace {
void set_env(const char* k, const char* v) {
#ifdef _WIN32
  _putenv_s(k, v ? v : "");
#else
  if (v) ::setenv(k, v, 1);
  else ::unsetenv(k);
#endif
}
const char* kManifest = R"({
 "format": "caps-release/1",
 "version": "0.1.0",
 "commit": "765DC9F3A1B2C3D4E5F60718293A4B5C6D7E8F90",
 "date": "2026-10-22T09:14:03Z",
 "title": "October 2026 update",
 "highlights": [
  {"title": "Semicrystalline stacks", "detail": "Lamellar cells with loops, bridges and tails.", "page": "semicrystalline"},
  {"title": "Thermal conductivity by MD"},
  "A plain line"
 ],
 "notes_url": "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest",
 "assets": {
  "macos-arm64": {"file": "CAPS-0.1.0-macos-arm64.dmg", "sha256": "AB12"},
  "windows-x64": {"file": "CAPS-0.1.0-windows-x64-setup.exe"},
  "linux-x86_64": "CAPS-0.1.0-linux-x86_64.AppImage"
 }
})";
ReleaseManifest release(const std::string& version, const std::string& commit, const std::string& date) {
  ReleaseManifest m;
  m.version = version, m.commit = commit, m.date = date;
  return m;
}
}  // namespace

TEST(Update, ParsesManifest) {
  const auto m = parse_release_manifest(kManifest);
  EXPECT_EQ(m.version, "0.1.0");
  EXPECT_EQ(m.commit, "765dc9f3a1b2c3d4e5f60718293a4b5c6d7e8f90");   // lower case
  EXPECT_EQ(m.title, "October 2026 update");
  ASSERT_EQ(m.highlights.size(), 3u);
  EXPECT_EQ(m.highlights[0].page, "semicrystalline");
  EXPECT_EQ(m.highlights[1].detail, "");
  EXPECT_EQ(m.highlights[2].title, "A plain line");
  ASSERT_EQ(m.assets.size(), 3u);
  EXPECT_EQ(m.assets.at("macos-arm64").sha256, "ab12");
  EXPECT_EQ(m.assets.at("linux-x86_64").file, "CAPS-0.1.0-linux-x86_64.AppImage");
  // written and read back the same
  const auto back = parse_release_manifest(release_manifest_json(m));
  EXPECT_EQ(back.commit, m.commit);
  EXPECT_EQ(back.highlights.size(), 3u);
  EXPECT_EQ(back.assets.at("macos-arm64").file, "CAPS-0.1.0-macos-arm64.dmg");
}

TEST(Update, RefusesOtherFiles) {
  EXPECT_THROW(parse_release_manifest("not json"), std::runtime_error);
  EXPECT_THROW(parse_release_manifest(R"({"tag_name": "v0.1.0"})"), std::runtime_error);
  EXPECT_THROW(parse_release_manifest(R"({"format": "caps-release/1", "version": "0.1.0", "date": "2026-10-22"})"), std::runtime_error);
}

TEST(Update, ComparesVersions) {
  EXPECT_EQ(compare_versions("0.1.0", "0.2.0"), -1);
  EXPECT_EQ(compare_versions("0.10.0", "0.9.9"), 1);
  EXPECT_EQ(compare_versions("v1.0", "1.0.0"), 0);
  EXPECT_EQ(compare_versions("0.2.0-rc1", "0.2.0"), -1);
  EXPECT_EQ(compare_versions("0.2.0+build5", "0.2.0"), 0);
}

TEST(Update, ReadsIsoDates) {
  EXPECT_EQ(parse_iso_time("1970-01-02"), 86400);
  EXPECT_EQ(parse_iso_time("2026-10-22T09:14:03Z"), parse_iso_time("2026-10-22T11:14:03+02:00"));
  EXPECT_EQ(parse_iso_time("2026-10-22T09:14:03.250Z"), parse_iso_time("2026-10-22T09:14:03Z"));
  EXPECT_LT(parse_iso_time("2026-10-21T23:00:00Z"), parse_iso_time("2026-10-22T00:30:00+01:00") + 1);
  EXPECT_EQ(parse_iso_time("unknown"), -1);
}

TEST(Update, SameCommitByPrefix) {
  EXPECT_TRUE(same_commit("765dc9f3a1", "765DC9F3A1B2C3D4E5F60718293A4B5C6D7E8F90"));
  EXPECT_FALSE(same_commit("765dc9f3a1", "96e9eb9a50"));
  EXPECT_FALSE(same_commit("unknown", "unknown"));
  EXPECT_FALSE(same_commit("765", "765dc9f3a1"));   // too short to tell
}

TEST(Update, HigherVersionIsAnUpdate) {
  const BuildIdentity local{"0.1.0", "96e9eb9a50", "2026-10-21T10:00:00Z"};
  // a higher version wins even with an older date; a lower one never is an update
  EXPECT_EQ(decide_update(local, release("0.2.0", "765dc9f3a1", "2026-01-01T00:00:00Z")), UpdateVerdict::Update);
  EXPECT_EQ(decide_update({"0.2.0", "96e9eb9a50", "2026-10-21T10:00:00Z"}, release("0.1.0", "765dc9f3a1", "2027-01-01")),
            UpdateVerdict::LocalNewer);
}

TEST(Update, EqualVersionsDecideByCommitAndDate) {
  const BuildIdentity local{"0.1.0", "96e9eb9a50", "2026-10-21T10:00:00+02:00"};
  // another commit, released later: an update
  EXPECT_EQ(decide_update(local, release("0.1.0", "765dc9f3a1b2c3d4", "2026-10-22T09:14:03Z")), UpdateVerdict::Update);
  // the same commit (short against full): none
  EXPECT_EQ(decide_update(local, release("0.1.0", "96e9eb9a50ffffffffffffffffffffffffffffff", "2026-10-22T09:14:03Z")), UpdateVerdict::Current);
  // this build is newer than the release: none
  EXPECT_EQ(decide_update(local, release("0.1.0", "765dc9f3a1", "2026-10-20T09:00:00Z")), UpdateVerdict::LocalNewer);
  // no date for this build: cannot tell
  EXPECT_EQ(decide_update({"0.1.0", "96e9eb9a50", "unknown"}, release("0.1.0", "765dc9f3a1", "2026-10-22")), UpdateVerdict::Unknown);
}

TEST(Update, EnvironmentTurnsItOff) {
  set_env("CAPS_NO_UPDATE_CHECK", nullptr);
  EXPECT_FALSE(update_checks_disabled());
  set_env("CAPS_NO_UPDATE_CHECK", "1");
  EXPECT_TRUE(update_checks_disabled());
  set_env("CAPS_NO_UPDATE_CHECK", "0");
  EXPECT_FALSE(update_checks_disabled());
  set_env("CAPS_NO_UPDATE_CHECK", "yes");
  EXPECT_TRUE(update_checks_disabled());
  set_env("CAPS_NO_UPDATE_CHECK", nullptr);
}

TEST(Update, LocalFeedAndDownloadUrl) {
  const auto dir = std::filesystem::temp_directory_path() / "caps_test_update";
  std::filesystem::create_directories(dir);
  const auto file = (dir / "release.json").string();
  std::ofstream(file) << kManifest;
  set_env("CAPS_UPDATE_FEED", file.c_str());
  EXPECT_EQ(update_feed(), file);
  const auto m = parse_release_manifest(fetch_release_manifest(update_feed()));
  set_env("CAPS_UPDATE_FEED", nullptr);
  EXPECT_EQ(update_feed(), kReleaseManifestUrl);
  EXPECT_EQ(download_url(m, "macos-arm64"), "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/CAPS-0.1.0-macos-arm64.dmg");
  EXPECT_EQ(download_url(m, "linux-x86_64-deb"), "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/CAPS-0.1.0-linux-x86_64.AppImage");
  EXPECT_EQ(download_url(m, "macos-x86_64"), "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest");
  EXPECT_THROW(fetch_release_manifest((dir / "missing.json").string()), std::runtime_error);
  std::filesystem::remove_all(dir);
}
