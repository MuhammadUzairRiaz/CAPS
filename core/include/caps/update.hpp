// Release manifests (caps-release/1) and the update decision: is a published CAPS newer than this build?
//
// A release carries release.json beside its installers; each app bundles the same file for its own build. The version
// decides first (a higher release version is an update, a lower one is not); with equal versions the build commit and
// its date decide (another commit with a later date is an update). Nothing here sends anything but the GET of the
// manifest itself: no cookies, no identifiers, no telemetry.
#pragma once

#include <map>
#include <string>
#include <vector>

namespace caps {

struct ReleaseHighlight {
  std::string title;
  std::string detail;
  std::string page;   // a Studio page to open ("semicrystalline", or a command id "module.open thermal"); may be empty
};

struct ReleaseAsset {
  std::string file;     // "CAPS-0.1.0-macos-arm64.dmg"
  std::string sha256;   // may be empty (the manifest bundled in an app is written before its installer exists)
};

struct ReleaseManifest {
  std::string version;   // "0.1.0"
  std::string commit;    // full SHA of the build commit
  std::string date;      // ISO 8601, the commit's date
  std::string title;
  std::vector<ReleaseHighlight> highlights;
  std::string notes_url;
  std::map<std::string, ReleaseAsset> assets;   // platform key ("macos-arm64", "windows-x64", "linux-x86_64-deb" …) → file
};

// Throws std::runtime_error on a file that is not caps-release/1 or lacks version, commit or date.
ReleaseManifest parse_release_manifest(const std::string& json);
std::string release_manifest_json(const ReleaseManifest& m);

// -1, 0, 1 for "0.1.0" < = > "0.2.0" (numeric parts; a pre-release "0.2.0-rc1" sorts before "0.2.0"); a missing part is 0.
int compare_versions(const std::string& a, const std::string& b);
// Seconds since 1970 (UTC) of "2026-10-22", "2026-10-22T09:14:03Z" or "…+02:00"; -1 when it cannot be read.
long long parse_iso_time(const std::string& s);
// Same commit: one is a prefix of the other (short and full SHAs), at least 7 characters, case ignored.
bool same_commit(const std::string& a, const std::string& b);

struct BuildIdentity {
  std::string version, commit, date;
};
BuildIdentity this_build();

enum class UpdateVerdict {
  Update,       // the release is newer: offer it
  Current,      // the release is this build
  LocalNewer,   // this build is newer than the release (a higher version, or a later commit)
  Unknown,      // equal versions, other commits, and a date that cannot be compared
};
UpdateVerdict decide_update(const BuildIdentity& local, const ReleaseManifest& remote);
const char* verdict_name(UpdateVerdict v);

// CAPS_NO_UPDATE_CHECK set to anything but "", "0", "false", "no" or "off": no update check, ever.
bool update_checks_disabled();
// The manifest of the latest release, and where it is fetched from (CAPS_UPDATE_FEED: a local file or another URL).
extern const char* const kReleaseManifestUrl;
std::string update_feed();
// The manifest's text: read from a file, or fetched with curl (HTTPS only, no config file, no cookies, a fixed
// user agent "CAPS", timeout_s seconds in all). Throws std::runtime_error when it cannot be had.
std::string fetch_release_manifest(const std::string& feed, int timeout_s = 6);
// This machine's platform key ("macos-arm64", "windows-x64", "linux-x86_64") and the download URL of the release's
// file for it (the release page when the manifest has none).
std::string platform_key();
std::string download_url(const ReleaseManifest& m, const std::string& key);

}  // namespace caps
